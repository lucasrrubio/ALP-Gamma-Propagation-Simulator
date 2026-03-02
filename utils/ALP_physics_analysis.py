import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
import os
import sys

# Noise-Free threshold for comparison
THRESHOLD_DIST = 800.0 

def generate_comparison_plots():
    if not os.path.exists("results/events.csv"):
        print("[ERROR] results/events.csv not found.")
        return

    df = pd.read_csv("results/events.csv")
    df['ebl_survival'] = np.exp(-df['tau_total'])
    df['boost'] = df['prob_survival'] / df['ebl_survival']
    
    # Calculate simulated fluxes for spectral fitting
    # Assuming intrinsic index 2.15 from your previous logs
    df['flux_sim'] = df['prob_survival'] * np.power(df['energy_tev'], -2.15)
    df['flux_theo'] = df['ebl_survival'] * np.power(df['energy_tev'], -2.15)

    fig, axs = plt.subplots(2, 2, figsize=(15, 12))
    plt.subplots_adjust(hspace=0.3, wspace=0.25)

    energy_bins = np.logspace(np.log10(0.1), np.log10(20.0), 15)

    for i, (is_filtered) in enumerate([False, True]):
        data = df[df['dist_mpc'] >= THRESHOLD_DIST].copy() if is_filtered else df.copy()
        label = f"Filtered (D > {THRESHOLD_DIST} Mpc)" if is_filtered else "Global (Including Nearby Sources)"
        
        # --- PREPARE STATISTICS ---
        # Group by energy and calculate mean and std
        grouped = data.groupby(pd.cut(data['energy_tev'], bins=energy_bins), observed=True)
        stats = grouped.agg({
            'boost': ['mean', 'std'],
            'flux_sim': ['mean', 'std'],
            'flux_theo': ['mean', 'std']
        }).dropna()
        
        x_centers = stats.index.map(lambda x: x.mid).values.astype(float)
        
        # --- TOP: BOOST FACTOR WITH STANDARD DEVIATION ---
        axs[0, i].errorbar(x_centers, stats['boost']['mean'], yerr=stats['boost']['std'], 
                           fmt='o', color='red', ecolor='red', elinewidth=1, capsize=2, 
                           alpha=0.7, label=r'Mean Boost $\pm 1\sigma$')
        
        axs[0, i].axhline(1.0, color='black', ls='--', alpha=0.6)
        axs[0, i].set_title(f"Boost Factor: {label}")
        axs[0, i].set_xscale('log')
        axs[0, i].set_yscale('log')
        axs[0, i].set_ylabel(r"Boost Factor ($F_{ALP}/F_{EBL}$)")
        axs[0, i].grid(True, which='both', alpha=0.2)
        axs[0, i].legend()

        # --- BOTTOM: SPECTRAL FIT WITH ERROR BARS AND SLOPE ERROR ---
        # Plot Theoretical EBL
        axs[1, i].errorbar(x_centers, stats['flux_theo']['mean'], yerr=stats['flux_theo']['std'],
                           fmt='s', color='black', ecolor='black', elinewidth=1, capsize=2, 
                           alpha=0.2, label='Standard EBL')
        
        # Plot Simulation with ALPs
        axs[1, i].errorbar(x_centers, stats['flux_sim']['mean'], yerr=stats['flux_sim']['std'],
                           fmt='o', color='blue', ecolor='blue', elinewidth=1, capsize=2, 
                           label='Simulation (ALPs)')

        # Fit with Uncertainty Estimation
        # Log-Log fit: log(Flux) = -Gamma * log(Energy) + Const
        log_x = np.log10(x_centers)
        
        # ALP Fit
        p_alp, cov_alp = np.polyfit(log_x, np.log10(stats['flux_sim']['mean']), 1, cov=True)
        gamma_alp = -p_alp[0]
        sigma_alp = np.sqrt(cov_alp[0,0])
        
        # EBL Fit
        p_ebl, cov_ebl = np.polyfit(log_x, np.log10(stats['flux_theo']['mean']), 1, cov=True)
        gamma_ebl = -p_ebl[0]
        sigma_ebl = np.sqrt(cov_ebl[0,0])
        
        # Calculate Hardening and Propagated Error
        delta_gamma = gamma_ebl - gamma_alp
        sigma_delta = np.sqrt(sigma_alp**2 + sigma_ebl**2)

        axs[1, i].set_title(rf"Spectral Fit: $\Delta\Gamma = {delta_gamma:.3f} \pm {sigma_delta:.3f}$")
        axs[1, i].set_xscale('log')
        axs[1, i].set_yscale('log')
        axs[1, i].set_xlabel("Energy [TeV]")
        axs[1, i].set_ylabel("Differential Flux")
        axs[1, i].grid(True, which='both', alpha=0.2)
        axs[1, i].legend()

    plt.suptitle("ALP Physics Analysis: Boost Factor and Spectral Hardening", fontsize=16)
    os.makedirs('results', exist_ok=True)
    plt.savefig('results/alp_comparison_scientific.png', dpi=300, bbox_inches='tight')
    print(f"[SUCCESS] Scientific comparison plot saved to results/")
    plt.show()

if __name__ == "__main__":
    generate_comparison_plots()