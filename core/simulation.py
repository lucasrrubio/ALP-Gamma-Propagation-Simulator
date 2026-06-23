import numpy as np
import pandas as pd
import multiprocessing as mp
from time import time
import json
from astropy.coordinates import SkyCoord
import astropy.units as u

from core.magnetic_fields import AGNJetField, IGMTurbulentField, GMFJanssonFarrar
from core.propagation_engine import (
    get_mixing_matrix,
    get_mixing_matrix_with_absorption,
    evolve_state,
    evolve_state_non_hermitian,
)
from utils.constants import (
    MPC_TO_INV_EV,
    GEV_INV_TO_INV_EV,
    KPC_TO_EV_INV,
    CM_TO_INV_EV,
)
from utils.ebl import EBLModel


def run_source_wrapper(args):
    source_data, config = args

    n_events = int(source_data.get("Number of Events", 1))

    m_a_sq = config["physics_parameters"]["m_a_eV"] ** 2
    g_ag   = config["physics_parameters"]["g_ag_GeV_inv"] * GEV_INV_TO_INV_EV

    ebl_model = EBLModel()
    igm_field = IGMTurbulentField(
        B_rms_nG=config["egmf_model"]["B_field_nG"],
        coherence_length_mpc=config["egmf_model"]["coherence_length_mpc"],
    )
    gmf_field = GMFJanssonFarrar()

    # Galactic coordinate handling
    if "GLAT" in source_data:
        l_coord, b_coord = source_data["GLON"], source_data["GLAT"]
    else:
        c = SkyCoord(
            ra=source_data.get("ra", 0) * u.degree,
            dec=source_data.get("dec", 0) * u.degree,
            frame="icrs",
        )
        l_coord, b_coord = c.galactic.l.degree, c.galactic.b.degree

    # Energy sampling (power law)
    gamma = config["source_generation"]["agn_energy_gamma"]
    E_min = config["source_generation"]["min_energy_tev"]
    E_max = config["source_generation"]["max_energy_tev"]
    r_vals      = np.random.random(n_events)
    energies_tev = (
        E_min ** (1 - gamma)
        + r_vals * (E_max ** (1 - gamma) - E_min ** (1 - gamma))
    ) ** (1 / (1 - gamma))

    # IGM step size
    dz_igm_inv_ev = config["egmf_model"]["coherence_length_mpc"] * MPC_TO_INV_EV

    # GMF step configuration
    dz_gmf_kpc    = 0.005
    dz_gmf_inv_ev = dz_gmf_kpc * KPC_TO_EV_INV

    # AGN field setup — parameters from HBL catalog (Tavecchio et al. 2010)
    B0_gauss = source_data.get("B0_G",    0.138)    # field at emission region [G]
    r0_cm    = source_data.get("r0_cm",   1.09e17)  # emission region distance [cm]
    delta_D  = source_data.get("delta_D", 25.0)     # Doppler factor of the jet

    agn_field = AGNJetField(B0_gauss=B0_gauss, r0_cm=r0_cm)

    # Log-spaced domains over [r0, 1000 r0], field sampled at geometric
    # midpoints. Resolution concentrates near r0 where B(r)=B0(r0/r) is
    # strongest and the oscillation fastest. Converges for all field strengths
    # and is ~7x faster than uniform linear stepping, which biases toward the
    # strong-field edge and fails to converge for B0 of a few G.
    N_AGN_STEPS   = 8000
    r_bounds_cm   = np.logspace(np.log10(r0_cm), np.log10(1000.0 * r0_cm),
                                N_AGN_STEPS)
    r_mid_cm      = np.sqrt(r_bounds_cm[:-1] * r_bounds_cm[1:])
    dz_agn_inv_ev = (r_bounds_cm[1:] - r_bounds_cm[:-1]) * CM_TO_INV_EV

    events = []
    for energy_tev in energies_tev:
        energy_ev = energy_tev * 1e12               # observer-frame energy
        energy_jet_ev = energy_ev / delta_D          # jet comoving-frame energy
        state     = np.array([1.0, 0.0, 0.0], dtype=np.complex128)

        # ── Region A: AGN jet ─────────────────────────────────────────────────
        # Mixing occurs in the jet comoving frame, where the B field is defined.
        # The photon energy is Doppler-shifted: E_jet = E_obs / delta_D.
        # The IGM and GMF blocks below use the observer-frame energy_ev.
        for k in range(len(r_mid_cm)):
            Bx, By = agn_field.get_field_vector(r_mid_cm[k])
            state  = evolve_state(
                state,
                get_mixing_matrix(energy_jet_ev, Bx, By, g_ag, m_a_sq),
                dz_agn_inv_ev[k],
            )

        # ── Region B: IGM with coherent EBL absorption ────────────────────────
        dist_mpc  = source_data.get("Distance", 100.0)
        tau_total = ebl_model.get_tau(energy_tev, dist_mpc)
        n_steps   = max(1, int(dist_mpc / config["egmf_model"]["coherence_length_mpc"]))

        gamma_per_step = (tau_total / n_steps) / dz_igm_inv_ev  # [eV]

        for _ in range(n_steps):
            Bx, By = igm_field.get_field_vector(new_cell=True)
            M_eff  = get_mixing_matrix_with_absorption(
                energy_ev, Bx, By, g_ag, m_a_sq, gamma_per_step
            )
            state = evolve_state_non_hermitian(state, M_eff, dz_igm_inv_ev)

        # ── Region C: Milky Way (GMF) ─────────────────────────────────────────
        d_kpc    = 20.0
        step_idx = 0
        Bt_e1, Bt_e2 = 0.0, 0.0   # stored turbulent components

        while d_kpc > 0:
            # Update turbulent cell every 20 steps (20 * 0.005 = 0.1 kpc)
            if step_idx % 20 == 0:
                Bt_e1, Bt_e2 = gmf_field.get_turbulent_field_transverse(
                    l_coord, b_coord, d_kpc
                )
            # Regular field recomputed at every step
            Br_e1, Br_e2 = gmf_field.get_regular_field_transverse(
                l_coord, b_coord, d_kpc
            )
            Bx, By = Br_e1 + Bt_e1, Br_e2 + Bt_e2
            state  = evolve_state(
                state,
                get_mixing_matrix(energy_ev, Bx, By, g_ag, m_a_sq),
                dz_gmf_inv_ev,
            )
            d_kpc    -= dz_gmf_kpc
            step_idx += 1

        events.append({
            "energy_tev":    float(energy_tev),
            "prob_survival": float(np.abs(state[0]) ** 2 + np.abs(state[1]) ** 2),
            "prob_alp":      float(np.abs(state[2]) ** 2),
            "dist_mpc":      float(dist_mpc),
            "tau_total":     float(tau_total),
        })

    return events


class SimulationManager:
    def __init__(self, config):
        self.config = config

    def run_full_scan(self, sources_df):
        n_src = len(sources_df)
        print(f"[PROCESS] Starting simulation for {n_src} sources...")
        start = time()

        tasks = [(src, self.config) for src in sources_df.to_dict("records")]
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = pool.map(run_source_wrapper, tasks)

        elapsed = time() - start
        avg     = elapsed / n_src if n_src > 0 else 0
        print(f"[TIME] Total execution time: {elapsed:.2f} s")
        print(f"[TIME] Average per source:   {avg:.4f} s")

        return pd.DataFrame([ev for sub in results for ev in sub])
