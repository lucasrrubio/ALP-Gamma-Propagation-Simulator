import numpy as np
import pandas as pd
import json
import os

def generate_catalog():
    with open("config.json", 'r') as f:
        config = json.load(f)
        
    n_src = config["simulation_parameters"]["n_sources"]
    print(f"[CATALOG] Generating {n_src} sources...")

    # Log-Normal distributions for B and Mass
    B_gauss = np.exp(np.random.normal(config["source_generation"]["agn_B_field_log_mean_G"], config["source_generation"]["agn_B_field_log_std_G"], n_src))
    M_sun = np.exp(np.random.normal(config["source_generation"]["agn_mass_log_mean_Msun"], config["source_generation"]["agn_mass_log_std_Msun"], n_src))
    Ledd = np.random.exponential(scale=config["source_generation"]["agn_ledd_exp_scale"], size=n_src)
    
    # Isotropic distribution on the sphere
    distance = np.random.uniform(config["simulation_parameters"]["min_distance_mpc"], config["simulation_parameters"]["max_distance_mpc"], n_src)
    glon = np.random.uniform(0, 360, n_src)
    glat = np.degrees(np.arcsin(np.random.uniform(-1, 1, n_src)))

    df = pd.DataFrame({
        'Bh': B_gauss, 'Mbh': M_sun, 'Ledd': Ledd,
        'Distance': distance, 'GLON': glon, 'GLAT': glat,
        'Number of Events': (config["simulation_parameters"]["total_events"] * Ledd // Ledd.sum()).astype(int)
    })
    
    output_path = config["file_paths"]["source_catalog_output"]
    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df[df['Number of Events'] > 0].to_csv(output_path, index=False)
    print(f"[SUCCESS] Catalog saved to {output_path}")

if __name__ == "__main__":
    generate_catalog()