import json
import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

def run_spectral_analysis(min_dist_mpc=800.0):
    print(f"--- Calculating Effective Spectral Index (D > {min_dist_mpc} Mpc) ---")
    
    if not os.path.exists("results/events.csv"):
        print("[ERROR] results/events.csv not found.")
        return
        
    df = pd.read_csv("results/events.csv")
    with open("config.json", 'r') as f:
        config = json.load(f)
    
    gamma_int = config["source_generation"]["agn_energy_gamma"]
    
    # Filter and prepare flux data (dN/dE ~ E^-gamma)
    df = df[df['dist_mpc'] >= min_dist_mpc].copy()
    df['flux_sim'] = df['prob_survival'] * np.power(df['energy_tev'], -gamma_int)
    df['flux_theo'] = np.exp(-df['tau_total']) * np.power(df['energy_tev'], -gamma_int)
    
    # Log-spaced binning for high-energy fit (1-20 TeV)
    bins = np.logspace(np.log10(1.0), np.log10(20.0), 15)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    df['ebin'] = pd.cut(df['energy_tev'], bins=bins, labels=bin_centers)
    stats = df.groupby('ebin', observed=True).agg({'flux_sim': 'mean', 'flux_theo': 'mean'}).dropna()
    
    x = stats.index.astype(float).values
    log_x = np.log10(x)
    
    # Linear regression in log-log space (log10(Flux) = -Gamma * log10(E) + C)
    coeff_sim = np.polyfit(log_x, np.log10(stats['flux_sim']), 1)
    gamma_eff_sim = -coeff_sim[0]
    
    coeff_theo = np.polyfit(log_x, np.log10(stats['flux_theo']), 1)
    gamma_eff_theo = -coeff_theo[0]

    # Visualization
    plt.figure(figsize=(10, 7))
    plt.loglog(x, stats['flux_sim'], 'ro', label='Simulation Data (ALPs)')
    plt.loglog(x, 10**coeff_sim[1] * x**coeff_sim[0], 'r-', 
               label=rf'ALP Fit: $\Gamma_{{obs}} = {gamma_eff_sim:.2f}$')
    
    plt.loglog(x, stats['flux_theo'], 'ks', label='Theory Data (EBL Only)', alpha=0.3)
    plt.loglog(x, 10**coeff_theo[1] * x**coeff_theo[0], 'k--', 
               label=rf'EBL Fit: $\Gamma_{{obs}} = {gamma_eff_theo:.2f}$')

    plt.xlabel('Energy [TeV]')
    plt.ylabel('Differential Flux [dN/dE]')
    plt.title(rf'Spectral Hardening Analysis ($D > {min_dist_mpc}$ Mpc)')
    plt.legend()
    plt.grid(True, which="both", alpha=0.2)
    plt.savefig('results/spectral_index_fit.png', dpi=300)
    
    print(f"\n--- PHYSICAL ANALYSIS RESULTS ---")
    print(f"Intrinsic Spectral Index: {gamma_int:.2f}")
    print(f"Observed Index (Standard EBL): {gamma_eff_theo:.2f} (Soft spectrum)")
    print(f"Observed Index (ALP Mixing): {gamma_eff_sim:.2f} (Hardened spectrum)")
    print(f"Spectral Hardening Signal (Delta Gamma): {gamma_eff_theo - gamma_eff_sim:.2f}")

if __name__ == "__main__":
    run_spectral_analysis()