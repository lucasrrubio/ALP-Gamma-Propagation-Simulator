import numpy as np
from utils.constants import NG_TO_EV2

class IGMTurbulentField:
    """
    Models the Intergalactic Magnetic Field (IGM) based on a turbulent cell approach.
    """
    def __init__(self, B_rms_nG, coherence_length_mpc):
        self.B_rms_nG = B_rms_nG
        self.coherence_length_mpc = coherence_length_mpc
        # Convert to natural units (eV^2)
        self.B_rms_ev2 = B_rms_nG * NG_TO_EV2

    def get_field_vector(self, new_cell=False):
        """
        Returns random transverse Bx, By components (Gaussian distribution).
        """
        if new_cell:
            # B_rms^2 = <Bx^2> + <By^2> + <Bz^2> -> sigma = B_rms / sqrt(3)
            sigma = self.B_rms_ev2 / np.sqrt(3.0)
            Bx = np.random.normal(0, sigma)
            By = np.random.normal(0, sigma)
            return Bx, By
        return 0.0, 0.0

class AGNJetField:
    """
    Models the AGN Jet Magnetic Field with radial power-law decay: B(r) = B_H * (R_H / r)^n.
    """
    def __init__(self, B_H_gauss, R_H_meter, n_index=2.0):
        # 1 Gauss = 1e9 nG
        self.B_H_ev2 = (B_H_gauss * 1e9) * NG_TO_EV2
        self.R_H = R_H_meter
        self.n = n_index
        # Coherent jet field angle
        self.psi = np.random.uniform(0, 2 * np.pi)

    def get_field_vector(self, r_meter):
        if r_meter < self.R_H:
            B_mag = self.B_H_ev2
        else:
            B_mag = self.B_H_ev2 * (self.R_H / r_meter)**self.n
            
        Bx = B_mag * np.cos(self.psi)
        By = B_mag * np.sin(self.psi)
        return Bx, By

class GMFJanssonFarrar:
    """
    Implements the Jansson-Farrar (2012) Galactic Magnetic Field (GMF) model.
    Includes Disk, Halo, and X-field regular components plus a turbulent disk component.
    """
    def __init__(self):
        # Regular components parameters (JF12)
        self.b_arms = np.array([0.1, 3.0, -0.9, -0.8, -2.0, -4.2, 0.0, 2.7]) 
        self.b_ring = 0.1 
        self.h_disk = 0.40 
        self.w_disk = 0.27 
        
        self.Bn = 1.4  
        self.Bs = -1.1 
        self.rh = 9.22 
        self.zh = 5.3  
        self.wh = 0.20 
        
        self.Bx0 = 4.6 
        self.ThetaX0 = np.radians(49.0)
        self.rc_x = 4.8 
        self.rx = 2.9   

        # Turbulent component parameters
        self.B_turb_rms = 0.1 # muG (RMS at center) 
        self.R_turb = 10.0    # kpc (scale radius)
        self.H_turb = 2.0     # kpc (scale height)

        self.R_sun = 8.5 
        self.muG_to_ev2 = 1000.0 * NG_TO_EV2

    def _logistic(self, x, h, w):
        return 1.0 / (1.0 + np.exp(-2.0 * (np.abs(x) - h) / w))

    def _get_disk_field(self, r, phi, z):
        if r > 20.0: return np.array([0.0, 0.0, 0.0])
        z_profile = (1.0 - self._logistic(z, self.h_disk, self.w_disk))
        B_cyl = 0.0
        if r < 3.0:
            B_cyl = 0.0 
        elif 3.0 <= r < 5.0:
            B_cyl = self.b_ring * (5.0 / r) * (1.0 - self._logistic(r, 3.0, 0.1))
        else:
            B_cyl = 2.0 * (5.0 / r) 
            
        return np.array([-B_cyl * np.sin(phi) * z_profile, B_cyl * np.cos(phi) * z_profile, 0.0])

    def _get_halo_field(self, r, z):
        z_profile = np.exp(-np.abs(z) / self.zh)
        r_profile = self._logistic(r, self.rh, self.wh)
        B_phi_val = self.Bn if z > 0 else self.Bs
        return B_phi_val * (1.0 - self._logistic(z, 0.0, self.wh)) * r_profile * z_profile

    def _get_x_field(self, r, phi, z):
        tan_theta = np.tan(self.ThetaX0)
        rp = max(0.1, r - np.abs(z) / tan_theta)
        B_x_mag = self.Bx0 * np.exp(-rp / self.rx) * (rp / r)
        sign_z = np.sign(z) if z != 0 else 1.0
        Br = B_x_mag * np.cos(self.ThetaX0) * sign_z 
        return np.array([Br * np.cos(phi), Br * np.sin(phi), B_x_mag * np.sin(self.ThetaX0)])

    def _get_turbulent_field(self, r, z):
        scale = np.exp(-r / self.R_turb) * np.exp(-np.abs(z) / self.H_turb)
        sigma = (self.B_turb_rms * scale) / np.sqrt(3.0)
        return np.random.normal(0, sigma, 3)

    def get_field_vector_transverse(self, l_deg, b_deg, d_kpc):
        l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)
        x_helio = d_kpc * np.cos(b_rad) * np.cos(l_rad)
        y_helio = d_kpc * np.cos(b_rad) * np.sin(l_rad)
        z_helio = d_kpc * np.sin(b_rad)
        
        x, y, z = x_helio - self.R_sun, y_helio, z_helio
        r, phi = np.sqrt(x**2 + y**2), np.arctan2(y, x)
        
        B_halo_mag = self._get_halo_field(r, z)
        B_halo = np.array([-B_halo_mag * np.sin(phi), B_halo_mag * np.cos(phi), 0.0])
        B_total = B_halo + self._get_x_field(r, phi, z) + self._get_disk_field(r, phi, z) + self._get_turbulent_field(r, z)
        
        e1 = np.array([-np.sin(b_rad) * np.cos(l_rad), -np.sin(b_rad) * np.sin(l_rad), np.cos(b_rad)])
        e2 = np.array([-np.sin(l_rad), np.cos(l_rad), 0.0])
        
        return np.dot(B_total, e1) * self.muG_to_ev2, np.dot(B_total, e2) * self.muG_to_ev2