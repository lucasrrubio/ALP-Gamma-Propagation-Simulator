import numpy as np
import pandas as pd
import multiprocessing as mp
from time import time
import json
from astropy.coordinates import SkyCoord
import astropy.units as u

from core.magnetic_fields import AGNJetField, IGMTurbulentField, GMFJanssonFarrar
from core.propagation_engine import get_mixing_matrix, evolve_state
from utils.constants import (
    MPC_TO_INV_EV, GEV_INV_TO_INV_EV, 
    G_CONST_SI, C_LIGHT_SI, M_SUN_KG
)
from utils.ebl import EBLModel

def run_source_wrapper(args):
    source_data, config = args
    source_id = source_data.get('source_id', 0)
    
    n_events_total = config["simulation_parameters"]["total_events"]
    n_sources = config["simulation_parameters"]["n_sources"]
    n_events = max(1, int(n_events_total / n_sources))

    m_a_sq = config["physics_parameters"]["m_a_eV"]**2
    g_ag = config["physics_parameters"]["g_ag_GeV_inv"] * GEV_INV_TO_INV_EV
    
    ebl_model = EBLModel() 
    igm_field = IGMTurbulentField(
        B_rms_nG=config["egmf_model"]["B_field_nG"],
        coherence_length_mpc=config["egmf_model"]["coherence_length_mpc"]
    )
    gmf_field = GMFJanssonFarrar()
    
    # Galactic Coordinate Handling
    if 'GLAT' in source_data:
        l_coord, b_coord = source_data['GLON'], source_data['GLAT']
    else:
        c = SkyCoord(ra=source_data.get('ra', 0)*u.degree, dec=source_data.get('dec', 0)*u.degree, frame='icrs')
        l_coord, b_coord = c.galactic.l.degree, c.galactic.b.degree

    # Energy Sampling (Power Law)
    gamma = config["source_generation"]["agn_energy_gamma"]
    E_min, E_max = config["source_generation"]["min_energy_tev"], config["source_generation"]["max_energy_tev"]
    r = np.random.random(n_events)
    energies_tev = (E_min**(1-gamma) + r*(E_max**(1-gamma) - E_min**(1-gamma)))**(1/(1-gamma))
    
    # Step configurations
    dz_agn_inv_ev = 1e-8 * MPC_TO_INV_EV
    dz_agn_meters = 1e-8 * 3.0857e22
    dz_igm_inv_ev = config["egmf_model"]["coherence_length_mpc"] * MPC_TO_INV_EV
    dz_gmf_inv_ev = 0.005 * 1e-3 * MPC_TO_INV_EV
    
    # AGN Field Setup
    rg = (G_CONST_SI * source_data.get('Mbh', 1e8) * M_SUN_KG) / (C_LIGHT_SI**2)
    agn_field = AGNJetField(B_H_gauss=source_data.get('Bh', 100.0), R_H_meter=rg)
    
    events = []
    for energy_tev in energies_tev:
        energy_ev = energy_tev * 1e12
        state = np.array([1.0, 0.0, 0.0], dtype=np.complex128)
        
        # Region A: AGN Jet
        r_curr, max_jet = 100.0 * rg, 1000.0 * rg
        while r_curr < max_jet:
            Bx, By = agn_field.get_field_vector(r_curr)
            state = evolve_state(state, get_mixing_matrix(energy_ev, Bx, By, g_ag, m_a_sq), dz_agn_inv_ev)
            r_curr += dz_agn_meters

        # Region B: IGM
        dist_mpc = source_data.get('Distance', 100.0)
        tau = ebl_model.get_tau(energy_tev, dist_mpc)
        n_steps = max(1, int(dist_mpc / config["egmf_model"]["coherence_length_mpc"]))
        att = np.exp(-(tau/n_steps)/2.0)
        for _ in range(n_steps):
            Bx, By = igm_field.get_field_vector(new_cell=True)
            state = evolve_state(state, get_mixing_matrix(energy_ev, Bx, By, g_ag, m_a_sq), dz_igm_inv_ev)
            state[0:2] *= att

        # Region C: Milky Way (GMF)
        d_kpc = 20.0
        while d_kpc > 0:
            Bx, By = gmf_field.get_field_vector_transverse(l_coord, b_coord, d_kpc)
            state = evolve_state(state, get_mixing_matrix(energy_ev, Bx, By, g_ag, m_a_sq), dz_gmf_inv_ev)
            d_kpc -= 0.005
            
        events.append({
            'source_id': int(source_id),
            'energy_tev': float(energy_tev),
            'prob_survival': float(np.abs(state[0])**2 + np.abs(state[1])**2),
            'prob_alp': float(np.abs(state[2])**2),
            'dist_mpc': float(dist_mpc),
            'tau_total': float(tau)
        })
    return events

class SimulationManager:
    def __init__(self, config):
        self.config = config
            
    def run_full_scan(self, sources_df):
        n_src = len(sources_df)
        print(f"[PROCESS] Starting simulation for {n_src} sources...")
        start = time()
        
        tasks = [(src, self.config) for src in sources_df.to_dict('records')]
        with mp.Pool(processes=mp.cpu_count()) as pool:
            results = pool.map(run_source_wrapper, tasks)
            
        elapsed = time() - start
        avg = elapsed / n_src if n_src > 0 else 0
        print(f"[TIME] Total execution time: {elapsed:.2f} s")
        print(f"[TIME] Average per source: {avg:.4f} s")
        
        return pd.DataFrame([ev for sub in results for ev in sub])