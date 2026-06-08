import numpy as np
from utils.constants import G_TO_EV2, NG_TO_EV2


class IGMTurbulentField:
    """
    Models the Intergalactic Magnetic Field (IGM) based on a turbulent cell approach.
    """
    def __init__(self, B_rms_nG, coherence_length_mpc):
        self.B_rms_nG = B_rms_nG
        self.coherence_length_mpc = coherence_length_mpc
        self.B_rms_ev2 = B_rms_nG * NG_TO_EV2

    def get_field_vector(self, new_cell=False):
        """
        Returns random transverse Bx, By components (Gaussian distribution).
        """
        if new_cell:
            sigma = self.B_rms_ev2 / np.sqrt(3.0)
            Bx = np.random.normal(0, sigma)
            By = np.random.normal(0, sigma)
            return Bx, By
        return 0.0, 0.0


class AGNJetField:
    """
    Models the AGN jet magnetic field with toroidal power-law decay.

    Field profile: B(r) = B0 * (r0 / r)   [toroidal, exponent n=1]

    The exponent n=1 reflects the dominance of the toroidal field component
    at distances r >= r0 from the jet base. The poloidal component decays as
    r^{-2} while the toroidal component decays as r^{-1}, so the toroidal
    term dominates beyond the emission region.

    Parameters
    ----------
    B0_gauss : float
        Magnetic field strength at the SSC emission region [G].
        Sampled from log-normal distribution calibrated on 14 HBLs
        from Tavecchio et al. (2010), MNRAS 401, 1570.
    r0_cm : float
        Distance of the emission region from the central black hole [cm],
        defined as r0 = R * delta_D from SSC one-zone fits.
        Sampled from log-normal distribution calibrated on the same sample.

    References
    ----------
    Gao et al. (2024), arXiv:2407.20118
    Pudritz et al. (2012), Space Sci. Rev. 169, 27
    Begelman, Blandford & Rees (1984), Rev. Mod. Phys. 56, 255
    """

    def __init__(self, B0_gauss, r0_cm):
        self.B0_ev2 = B0_gauss * G_TO_EV2   # convert G to eV^2
        self.r0     = r0_cm                  # [cm]
        self.psi    = np.random.uniform(0, 2 * np.pi)  # random field orientation

    def get_field_vector(self, r_cm):
        """
        Returns (Bx, By) transverse field components in eV^2 at distance r_cm.
        Clamps to r0 to avoid extrapolation inside the emission region.
        """
        r     = max(r_cm, self.r0)
        B_mag = self.B0_ev2 * (self.r0 / r)   # toroidal decay: B proportional to r^{-1}
        return B_mag * np.cos(self.psi), B_mag * np.sin(self.psi)


class GMFJanssonFarrar:
    """
    Implements the Jansson-Farrar (2012) Galactic Magnetic Field (GMF) model.
    Includes Disk, Halo, and X-field regular components plus a turbulent disk component.
    """
    def __init__(self):
        # Regular disk parameters (JF12)
        self.b_arms = np.array([0.1, 3.0, -0.9, -0.8, -2.0, -4.2, 0.0, 2.7])
        self.r_arms = np.array([5.1, 6.3, 7.1, 8.3, 9.8, 11.4, 12.7, 15.5])  # kpc
        self.pitch_angle = np.radians(11.5)
        self.b_ring = 0.1
        self.h_disk = 0.40
        self.w_disk = 0.27

        # Halo parameters
        self.Bn = 1.4
        self.Bs = -1.1
        self.rh = 9.22
        self.zh = 5.3
        self.wh = 0.20

        # X-field parameters
        self.Bx0     = 4.6
        self.ThetaX0 = np.radians(49.0)
        self.rc_x    = 4.8
        self.rx      = 2.9

        # Turbulent component parameters
        self.B_turb_rms = 0.1   # muG (RMS at centre)
        self.R_turb     = 10.0  # kpc (scale radius)
        self.H_turb     = 2.0   # kpc (scale height)

        # Global constants
        self.R_sun      = 8.5
        self.muG_to_ev2 = 1000.0 * NG_TO_EV2

    def _logistic(self, x, h, w):
        return 1.0 / (1.0 + np.exp(-2.0 * (np.abs(x) - h) / w))

    def _get_disk_field(self, r, phi, z):
        """Calculates the regular disk component in muG."""
        if r > 20.0 or r < 3.0:
            return np.array([0.0, 0.0, 0.0])
        z_profile = (1.0 - self._logistic(z, self.h_disk, self.w_disk))
        if 3.0 <= r < 5.0:
            B_cyl = self.b_ring * (5.0 / r) * (1.0 - self._logistic(r, 3.0, 0.1))
        else:
            inv_tan_i = 1.0 / np.tan(self.pitch_angle)
            phase     = (phi - inv_tan_i * np.log(r / self.R_sun)) % (2 * np.pi)
            idx       = int((phase / (2 * np.pi)) * 8) % 8
            B_cyl     = self.b_arms[idx] * (5.0 / r)
        return np.array([-B_cyl * np.sin(phi) * z_profile,
                          B_cyl * np.cos(phi) * z_profile,
                          0.0])

    def _get_halo_field(self, r, z):
        """Calculates the toroidal halo component magnitude in muG."""
        halo_vertical_transition = self._logistic(z, self.h_disk, self.w_disk)
        z_profile   = np.exp(-np.abs(z) / self.zh)
        r_profile   = 1.0 - self._logistic(r, self.rh, self.wh)
        B_phi_val   = self.Bn if z > 0 else self.Bs
        return B_phi_val * halo_vertical_transition * r_profile * z_profile

    def _get_x_field(self, r, phi, z):
        """Calculates the poloidal X-field component in muG."""
        if r == 0:
            return np.array([0.0, 0.0, 0.0])
        tan_theta0 = np.tan(self.ThetaX0)
        rp = r - np.abs(z) / tan_theta0
        if rp < 0:
            return np.array([0.0, 0.0, 0.0])
        if rp < self.rc_x:
            theta_x = np.arctan(tan_theta0 * rp / self.rc_x)
            B_x_mag = self.Bx0 * np.exp(-rp / self.rx) * (self.rc_x / r) ** 2
        else:
            theta_x = self.ThetaX0
            B_x_mag = self.Bx0 * np.exp(-rp / self.rx) * (rp / r)
        sign_z = np.sign(z) if z != 0 else 1.0
        Br = B_x_mag * np.cos(theta_x) * sign_z
        Bz = B_x_mag * np.sin(theta_x)
        return np.array([Br * np.cos(phi), Br * np.sin(phi), Bz])

    def _get_turbulent_field(self, r, z):
        """Calculates the turbulent disk component in muG."""
        scale = np.exp(-r / self.R_turb) * np.exp(-np.abs(z) / self.H_turb)
        sigma = (self.B_turb_rms * scale) / np.sqrt(3.0)
        return np.random.normal(0, sigma, 3)

    def get_field_vector_transverse(self, l_deg, b_deg, d_kpc):
        """
        Transforms coordinates and returns transverse B components in eV^2.
        """
        l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)
        x_helio = d_kpc * np.cos(b_rad) * np.cos(l_rad)
        y_helio = d_kpc * np.cos(b_rad) * np.sin(l_rad)
        z_helio = d_kpc * np.sin(b_rad)
        x, y, z = x_helio - self.R_sun, y_helio, z_helio
        r, phi  = np.sqrt(x**2 + y**2), np.arctan2(y, x)
        B_halo_mag = self._get_halo_field(r, z)
        B_halo     = np.array([-B_halo_mag * np.sin(phi),
                                B_halo_mag * np.cos(phi),
                                0.0])
        B_total = (B_halo
                   + self._get_x_field(r, phi, z)
                   + self._get_disk_field(r, phi, z)
                   + self._get_turbulent_field(r, z))
        e1 = np.array([-np.sin(b_rad) * np.cos(l_rad),
                       -np.sin(b_rad) * np.sin(l_rad),
                        np.cos(b_rad)])
        e2 = np.array([-np.sin(l_rad), np.cos(l_rad), 0.0])
        return (np.dot(B_total, e1) * self.muG_to_ev2,
                np.dot(B_total, e2) * self.muG_to_ev2)

    def get_regular_field_transverse(self, l_deg, b_deg, d_kpc):
        """Returns only the regular GMF components (disk + halo + X-field) in eV^2."""
        l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)
        x_helio = d_kpc * np.cos(b_rad) * np.cos(l_rad)
        y_helio = d_kpc * np.cos(b_rad) * np.sin(l_rad)
        z_helio = d_kpc * np.sin(b_rad)
        x, y, z = x_helio - self.R_sun, y_helio, z_helio
        r, phi  = np.sqrt(x**2 + y**2), np.arctan2(y, x)
        B_halo_mag = self._get_halo_field(r, z)
        B_halo     = np.array([-B_halo_mag * np.sin(phi),
                                B_halo_mag * np.cos(phi),
                                0.0])
        B_reg = (B_halo
                 + self._get_x_field(r, phi, z)
                 + self._get_disk_field(r, phi, z))
        e1 = np.array([-np.sin(b_rad) * np.cos(l_rad),
                       -np.sin(b_rad) * np.sin(l_rad),
                        np.cos(b_rad)])
        e2 = np.array([-np.sin(l_rad), np.cos(l_rad), 0.0])
        return (np.dot(B_reg, e1) * self.muG_to_ev2,
                np.dot(B_reg, e2) * self.muG_to_ev2)

    def get_turbulent_field_transverse(self, l_deg, b_deg, d_kpc):
        """Returns only the turbulent GMF component in eV^2."""
        l_rad, b_rad = np.radians(l_deg), np.radians(b_deg)
        x_helio = d_kpc * np.cos(b_rad) * np.cos(l_rad)
        y_helio = d_kpc * np.cos(b_rad) * np.sin(l_rad)
        z_helio = d_kpc * np.sin(b_rad)
        x, y, z = x_helio - self.R_sun, y_helio, z_helio
        r, phi  = np.sqrt(x**2 + y**2), np.arctan2(y, x)
        B_turb  = self._get_turbulent_field(r, z)
        e1 = np.array([-np.sin(b_rad) * np.cos(l_rad),
                       -np.sin(b_rad) * np.sin(l_rad),
                        np.cos(b_rad)])
        e2 = np.array([-np.sin(l_rad), np.cos(l_rad), 0.0])
        return (np.dot(B_turb, e1) * self.muG_to_ev2,
                np.dot(B_turb, e2) * self.muG_to_ev2)