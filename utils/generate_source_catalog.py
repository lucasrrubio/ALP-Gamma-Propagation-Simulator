import numpy as np
import pandas as pd
import json
import os

# ── Log-normal parameters calibrated on 14 standard HBLs ─────────────────────
# Source: Tavecchio et al. (2010), MNRAS 401, 1570
# Extreme-HBLs (B0 < 0.01 G) excluded: 1ES 0229+200, 1ES 0347-121, B3 2322+396
# Fermi-LAT fluxes: Wang et al. (2023), arXiv:2312.01122
# All mu/sigma values in natural logarithm (ln) base

# B0 [G] — magnetic field at the SSC emission region
B0_MU_LN  = -1.9824
B0_SIG_LN =  0.9938

# R [cm] — radius of the SSC emission region
R_MU_LN  = 36.1369
R_SIG_LN =  0.6372

# delta_D — Doppler factor of the emitting plasma
# Used both to derive r0 = R * delta_D and for the relativistic Doppler
# boost of the photon energy in the jet frame (E_jet = E_obs / delta_D).
# Note: delta_D takes discrete values in Tavecchio et al. (2010); the
# continuous log-normal is an idealisation for smooth population sampling.
DELTA_MU_LN  = 3.2947
DELTA_SIG_LN = 0.2285

# F1000 [ph/cm^2/s] — Fermi-LAT integrated flux 1-100 GeV (weighting only)
FLUX_MU_LN  = -19.2716
FLUX_SIG_LN =   1.3481


def generate_catalog():
    with open("config.json", "r") as f:
        config = json.load(f)

    n_src = config["simulation_parameters"]["n_sources"]
    print(f"[CATALOG] Generating {n_src} synthetic HBL sources...")
    print(f"[CATALOG] Parameters calibrated on Tavecchio et al. (2010), N=14")
    print(f"[CATALOG] r0 = R * delta_D (Form A); Doppler boost E_jet = E_obs/delta_D")
    print(f"[CATALOG] Weighting: Fermi-LAT Flux1000 (Wang et al. 2023)")

    rng = np.random.default_rng()

    # ── Sample B0 from calibrated log-normal distribution ────────────────────
    B0_G  = np.exp(rng.normal(B0_MU_LN, B0_SIG_LN, n_src))   # [G]

    # ── Sample R and delta_D independently, then derive r0 = R * delta_D ─────
    # Pearson correlation between log10(R) and log10(delta_D) in the calibration
    # sample is -0.021 (negligible), justifying independent sampling.
    R_cm    = np.exp(rng.normal(R_MU_LN,     R_SIG_LN,     n_src))   # [cm]
    delta_D = np.exp(rng.normal(DELTA_MU_LN, DELTA_SIG_LN, n_src))   # dimensionless
    r0_cm   = R_cm * delta_D                                          # [cm]

    # ── Sample F1000 from calibrated log-normal distribution ─────────────────
    # Used exclusively as a luminosity weight — does not enter propagation
    flux1000 = np.exp(rng.normal(FLUX_MU_LN, FLUX_SIG_LN, n_src))  # [ph/cm^2/s]

    # ── Isotropic sky distribution ────────────────────────────────────────────
    distance = rng.uniform(
        config["simulation_parameters"]["min_distance_mpc"],
        config["simulation_parameters"]["max_distance_mpc"],
        n_src
    )
    glon = rng.uniform(0, 360, n_src)
    glat = np.degrees(np.arcsin(rng.uniform(-1, 1, n_src)))

    # ── Event weighting proportional to F1000 ────────────────────────────────
    target_events = config["simulation_parameters"]["total_events"]
    n_events_raw  = (flux1000 / flux1000.sum() * target_events).astype(int)
    n_events      = np.clip(n_events_raw, 1, None)

    df = pd.DataFrame({
        "B0_G":             B0_G,
        "R_cm":             R_cm,
        "delta_D":          delta_D,
        "r0_cm":            r0_cm,
        "flux1000":         flux1000,
        "Distance":         distance,
        "GLON":             glon,
        "GLAT":             glat,
        "Number of Events": n_events,
    })

    output_path = config["file_paths"]["source_catalog_output"]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df.to_csv(output_path, index=False)

    print(f"[SUCCESS] Catalog saved to {output_path}")
    print(f"[INFO]    B0:      median={np.median(B0_G):.3f} G, "
          f"range=[{B0_G.min():.3e}, {B0_G.max():.3e}] G")
    print(f"[INFO]    R:       median={np.median(R_cm):.3e} cm, "
          f"range=[{R_cm.min():.3e}, {R_cm.max():.3e}] cm")
    print(f"[INFO]    delta_D: median={np.median(delta_D):.1f}, "
          f"range=[{delta_D.min():.1f}, {delta_D.max():.1f}]")
    print(f"[INFO]    r0:      median={np.median(r0_cm):.3e} cm, "
          f"range=[{r0_cm.min():.3e}, {r0_cm.max():.3e}] cm")
    print(f"[INFO]    Events:  total={n_events.sum()}, "
          f"min={n_events.min()}, max={n_events.max()}")


if __name__ == "__main__":
    generate_catalog()
