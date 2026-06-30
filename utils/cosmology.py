"""
Cosmology helpers for the redshift-resolved IGM treatment.

The intergalactic-medium propagation can be run either as a uniform slab at
z = 0 (legacy) or with the background evolved along the line of sight. The
cosmological grid below reproduces gammaALPs' ``trafo.cosmo_cohlength`` exactly
(cell grid agrees to ~1e-6 in z and dL for z = 0.433): cells of observer-frame
coherence length dL = L0 / (1 + z') spaced by luminosity distance in a flat
LambdaCDM cosmology (H0 = 70 km/s/Mpc, Om = 0.3, OL = 0.7).

This module is self-contained (no dependency on the validation apparatus) so the
production engine can import it directly.
"""

import numpy as np

# Flat LambdaCDM background (matches gammaALPs' default FlatLambdaCDM)
H0_KM_S_MPC = 70.0
C_KM_S      = 299792.458
OMEGA_M     = 0.3
OMEGA_L     = 0.7

# Hubble distance c/H0 [Mpc]; also the linear redshift<->distance scale that the
# tabulated EBL optical depth tau(z, E) assumes (z = distance_mpc * H0 / c).
HUBBLE_DIST_MPC = C_KM_S / H0_KM_S_MPC


def redshift_from_distance(distance_mpc):
    """Redshift implied by a Hubble-flow distance, the same convention the EBL
    table uses internally (z = distance * H0 / c)."""
    return distance_mpc * H0_KM_S_MPC / C_KM_S


def cosmo_cohlength(zmax, L0_mpc, Om=OMEGA_M, OL=OMEGA_L):
    """IGM coherence-length / redshift grid, replicating gammaALPs'
    trafo.cosmo_cohlength.

    Cells of observer-frame length dL = L0 / (1 + z') are placed by luminosity
    distance in a flat LambdaCDM cosmology. The last cell is trimmed so the grid
    terminates exactly at the source.

    Parameters
    ----------
    zmax : float
        Source redshift.
    L0_mpc : float
        Observer-frame coherence length [Mpc].

    Returns
    -------
    (dL_mpc, z_step) : tuple of np.ndarray
        dL_mpc[N]   : observer-frame cell lengths [Mpc].
        z_step[N+1] : cell-edge redshifts, ascending from 0 (observer) to zmax
                      (source). The per-cell local redshift is z_step[:-1] (the
                      lower edge of each cell), matching gammaALPs.
    """
    zg = np.linspace(0.0, zmax * 1.1, 4000)
    Ez = np.sqrt(Om * (1.0 + zg) ** 3 + OL)
    Dc = HUBBLE_DIST_MPC * np.concatenate(
        [[0.0], np.cumsum(0.5 * (1.0 / Ez[1:] + 1.0 / Ez[:-1]) * np.diff(zg))])
    Dlum = (1.0 + zg) * Dc                      # luminosity distance [Mpc]
    zfunc = lambda lum: float(np.interp(lum, Dlum, zg))

    z = [0.0, zfunc(L0_mpc)]
    dL = [L0_mpc]
    while z[-1] < zmax:
        dL.append(L0_mpc / (1.0 + z[-1]))
        z.append(zfunc(sum(dL)))
    dL[-1] = float(np.interp(zmax, zg, Dlum) - np.interp(z[-2], zg, Dlum))
    z[-1] = zmax
    return np.array(dL), np.array(z)
