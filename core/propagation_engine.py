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

_unitarity_violation_count = 0

def get_mixing_matrix_with_absorption(energy_ev, Bx_ev2, By_ev2,
                                       g_ag_inv_ev, m_a_sq_ev2,
                                       gamma_absorb_ev):
    delta_agx = 0.5 * g_ag_inv_ev * Bx_ev2
    delta_agy = 0.5 * g_ag_inv_ev * By_ev2
    delta_a   = m_a_sq_ev2 / (2.0 * energy_ev)

    abs_diag = -0.5j * gamma_absorb_ev

    return np.array([
        [abs_diag,    0.0,         delta_agx],
        [0.0,         abs_diag,    delta_agy],
        [delta_agx,   delta_agy,   delta_a  ]
    ], dtype=np.complex128)


def evolve_state(psi, M, dist_step_ev_inv):
    """
    Evolves the state vector psi via matrix exponentiation.
    M is real symmetric (Hermitian), so eigh is exact and returns
    real eigenvalues — no clipping needed, no phase loss.
    """
    global _unitarity_violation_count
    try:
        norm_in = np.real(np.dot(psi.conj(), psi))

        w, v = np.linalg.eigh(M)
        phase = np.exp(-1j * w * dist_step_ev_inv)
        psi_out = v @ (phase * (v.conj().T @ psi))

        # Unitarity validation
        norm_out = np.real(np.dot(psi_out.conj(), psi_out))
        if norm_in > 1e-30 and abs(norm_out - norm_in) / norm_in > 1e-8:
            _unitarity_violation_count += 1
            if _unitarity_violation_count <= 10:
                import sys
                print(
                    f"[WARN] Unitarity violation in step: "
                    f"norm_in²={norm_in:.10f}  norm_out²={norm_out:.10f}  "
                    f"delta={abs(norm_out-norm_in):.2e}  "
                    f"(total violations: {_unitarity_violation_count})",
                    file=sys.stderr
                )
 
        return psi_out

    except np.linalg.LinAlgError as e:
        import sys
        print(f"[ERROR] eigh failed: {e} — returning unnormalized psi", file=sys.stderr)
        return psi

def evolve_state_non_hermitian(psi, M, dist_step_ev_inv):
    try:
        A = -1j * M * dist_step_ev_inv
        w, v = np.linalg.eig(A)
        exp_w = np.exp(w)
        b = np.linalg.solve(v, psi)
        psi_out = v @ (exp_w * b)

        norm_sq = np.real(np.dot(psi_out.conj(), psi_out))
        if norm_sq > 1.0 + 1e-9:
            import sys
            print(f"[ERROR] norm²={norm_sq:.6f} > 1 — negative absorption, verify EBL gamma",
                  file=sys.stderr)
        return psi_out
    except np.linalg.LinAlgError:
        return psi