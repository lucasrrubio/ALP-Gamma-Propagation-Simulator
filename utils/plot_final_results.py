import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os

# Scientific Plotting Configuration
plt.style.use('seaborn-v0_8-paper')
plt.rcParams.update({
    "font.family": "serif",
    "axes.labelsize": 14,
    "font.size": 12,
    "legend.fontsize": 12,
    "xtick.labelsize": 12,
    "ytick.labelsize": 12,
    "figure.figsize": (10, 7)
})

def plot_filtered_results(file_path="results/events.csv", min_dist_mpc=800.0):
    print(f"--- Generating Analytical Plots (Distance > {min_dist_mpc} Mpc) ---")
    
    if not os.path.exists(file_path):
        print(f"[ERROR] Result file {file_path} not found.")
        return
        
    df = pd.read_csv(file_path)
    
    # Filter for distant sources to isolate ALP-induced transparency
    mask = df['dist_mpc'] >= min_dist_mpc
    df_filtered = df[mask].copy()
    
    if len(df_filtered) == 0:
        print("[ERROR] No sources found above the distance threshold.")
        return
        
    print(f"[INFO] Processing {len(df_filtered)} filtered events for ALP physics analysis.")

    # Calculate theoretical survival (EBL only)
    df_filtered['prob_theoretical'] = np.exp(-df_filtered['tau_total'])

    # Logarithmic energy binning (0.1 to 20 TeV)
    bins = np.logspace(np.log10(0.1), np.log10(20.0), 25)
    bin_centers = (bins[:-1] + bins[1:]) / 2
    df_filtered['energy_bin'] = pd.cut(df_filtered['energy_tev'], bins=bins, labels=bin_centers)
    
    # Statistical aggregation
    grouped = df_filtered.groupby('energy_bin', observed=True)
    mean_sim = grouped['prob_survival'].mean()
    mean_theo = grouped['prob_theoretical'].mean()
    
    valid = ~mean_sim.isna()
    x = mean_sim.index.astype(float)[valid]
    y_sim = mean_sim[valid]
    y_theo = mean_theo[valid]
    
    # --- Figure 1: Survival Probability (ALP Signal vs Standard Physics) ---
    plt.figure()
    plt.plot(x, y_theo, 'k--', linewidth=2, label='Standard Physics (EBL Only)', alpha=0.7)
    plt.plot(x, y_sim, 'r-o', linewidth=2, label='ALP-Photon Mixing (Simulation)', markersize=6)
    
    plt.yscale('log')
    plt.xscale('log')
    plt.ylim(1e-4, 1.1)
    plt.ylabel(r'Survival Probability ($P_{\gamma\gamma}$)')
    plt.xlabel('Energy [TeV]')
    plt.title(f'Photon Survival Probability ($D > {min_dist_mpc}$ Mpc)')
    plt.grid(True, which="both", alpha=0.2)
    plt.legend()
    plt.savefig('results/final_spectrum_filtered.png', dpi=300)
    print("-> Saved: results/final_spectrum_filtered.png")
    plt.close()
    
    # --- Figure 2: Boost Factor (Anomalous Transparency Signal) ---
    boost = y_sim / np.maximum(y_theo, 1e-10) 
    
    plt.figure()
    plt.plot(x, boost, 'b-o', linewidth=2, label='Boost Factor')
    plt.axhline(1.0, color='r', linestyle='--', label='Standard Model Reference')
    
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Boost Factor ($P_{ALP} / P_{EBL}$)')
    plt.xlabel('Energy [TeV]')
    plt.title(f'ALP-induced Flux Amplification ($D > {min_dist_mpc}$ Mpc)')
    plt.grid(True, which="both", alpha=0.2)
    plt.legend()
    plt.savefig('results/final_boost_filtered.png', dpi=300)
    print("-> Saved: results/final_boost_filtered.png")
    
    # Identify Maximum Boost
    high_e_mask = x > 8.0 
    if any(high_e_mask):
        max_b = boost[high_e_mask].max()
        print(f"\n[RESULT] MAXIMUM BOOST FACTOR DETECTED: {max_b:.2f}x")
    plt.close()

if __name__ == "__main__":
    plot_filtered_results(min_dist_mpc=800.0)