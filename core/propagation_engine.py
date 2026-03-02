import numpy as np

def get_mixing_matrix(energy_ev, Bx_ev2, By_ev2, g_ag_inv_ev, m_a_sq_ev2):
    """
    Constructs the 3x3 photon-ALP mixing matrix (M).
    Units: [eV]
    """
    delta_agx = 0.5 * g_ag_inv_ev * Bx_ev2
    delta_agy = 0.5 * g_ag_inv_ev * By_ev2
    delta_a = m_a_sq_ev2 / (2.0 * energy_ev)
    
    return np.array([
        [0.0,         0.0,         delta_agx],
        [0.0,         0.0,         delta_agy],
        [delta_agx,   delta_agy,   delta_a]
    ], dtype=np.complex128)

def evolve_state(psi, M, dist_step_ev_inv):
    """
    Evolves the state vector psi using stable diagonalization and Safe Clipping for unitarity.
    """
    try:
        A = -1j * M * dist_step_ev_inv
        w, v = np.linalg.eig(A)
        # Numerical stability: clip real part to avoid probability explosion
        exp_w = np.exp(np.clip(np.real(w), -700.0, 0.0)) * np.exp(1j * np.imag(w))
        b = np.linalg.solve(v, psi)
        return v @ (exp_w * b)
    except (np.linalg.LinAlgError, ValueError):
        return psi