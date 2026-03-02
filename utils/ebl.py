import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator

class EBLModel:
    """
    Realistic EBL model implementation (Dominguez 2011).
    Parses .out grid files to interpolate optical depth (tau).
    """
    def __init__(self):
        base_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.file_path = os.path.join(base_dir, "data", "ebl_models", "tau_dominguez11_cta.out")
        
        if not os.path.exists(self.file_path):
            raise FileNotFoundError(f"EBL data not found at: {self.file_path}")
        self._load_grid()

    def _load_grid(self):
        full_data = np.loadtxt(self.file_path, comments='#')
        with open(self.file_path, 'r') as f:
            for line in f:
                if line.strip() and not line.startswith('#'):
                    energies_tev = np.fromstring(line.strip(), sep=' ')
                    if energies_tev[0] == 0.0: energies_tev = energies_tev[1:]
                    break
        
        matrix_data = full_data[1:] if np.allclose(full_data[0, -len(energies_tev):], energies_tev) else full_data
        self.redshifts = matrix_data[:, 0]
        self.interpolator = RegularGridInterpolator(
            (self.redshifts, np.log10(energies_tev)), matrix_data[:, 1:], 
            bounds_error=False, fill_value=None 
        )

    def get_tau(self, energy_tev, distance_mpc):
        z = max(1e-5, distance_mpc * 70.0 / 299792.458)
        if energy_tev <= 0: return 0.0
        tau = self.interpolator(np.array([z, np.log10(energy_tev)]))
        return float(max(0.0, tau.item() if isinstance(tau, np.ndarray) else tau))