import os
import numpy as np
from scipy.interpolate import RegularGridInterpolator


class EBLModel:
    """
    EBL opacity model: Dominguez et al. (2011), MNRAS 410, 2556.

    Data source: ebltable package (tau_from_model.OptDepth),
    exported to tau_dominguez11_full.out.

    This file replaces the previous tau_dominguez11_cta.out which produced
    incorrect (underestimated) tau values at E > 3 TeV due to the limited
    energy grid of the CTA-optimised version. The new file covers
    1 GeV -- 100 TeV and z = 0.001 -- 10, matching the gammaALPs backend.

    File format:
        Row 0 : [0, log10(E1/TeV), log10(E2/TeV), ...]
        Rows 1+: [z, tau(E1), tau(E2), ...]
    """

    def __init__(self):
        base_dir       = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
        self.file_path = os.path.join(
            base_dir, "data", "ebl_models", "tau_dominguez11_full.out"
        )

        if not os.path.exists(self.file_path):
            raise FileNotFoundError(
                f"EBL data not found at: {self.file_path}\n"
                f"Run find_ebl_file.py to generate it from the ebltable package."
            )
        self._load_grid()

    def _load_grid(self):
        data = np.loadtxt(self.file_path, comments='#')

        # Row 0: [0, log10(E1), log10(E2), ...]
        energies_log10 = data[0, 1:]     # log10(E/TeV)

        # Rows 1+: [z, tau(E1), tau(E2), ...]
        self.redshifts  = data[1:, 0]
        tau_matrix      = data[1:, 1:]

        self.interpolator = RegularGridInterpolator(
            (self.redshifts, energies_log10),
            tau_matrix,
            bounds_error=False,
            fill_value=None,
        )

    def get_tau(self, energy_tev, distance_mpc):
        """
        Returns the EBL optical depth tau(E, z).

        Parameters
        ----------
        energy_tev : float
            Photon energy in TeV.
        distance_mpc : float
            Source distance in Mpc (converted to redshift via Hubble flow).

        Returns
        -------
        float
            Optical depth tau >= 0.
        """
        z = max(1e-5, distance_mpc * 70.0 / 299792.458)
        if energy_tev <= 0:
            return 0.0
        tau = self.interpolator(np.array([z, np.log10(energy_tev)]))
        return float(max(0.0, tau.item() if isinstance(tau, np.ndarray) else tau))
