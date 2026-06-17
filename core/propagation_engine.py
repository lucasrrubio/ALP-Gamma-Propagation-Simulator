import numpy as np

# ── QED vacuum birefringence constants ────────────────────────────────────────
# The QED vacuum birefringence (Heisenberg-Euler effective Lagrangian) modifies
# the photon dispersion in a magnetic field. Crucially, the two photon
# polarization states acquire DIFFERENT refractive indices:
#   - parallel to B:      Delta_par_QED  = 3.5 * DQED   (coefficient 7/2)
#   - perpendicular to B: Delta_perp_QED = 2.0 * DQED   (coefficient 4/2)
# Only the photon component PARALLEL to the transverse B field couples to the
# ALP. This is implemented by working in the B-aligned basis and rotating to
# the lab frame by the field angle psi.
#
# References:
#   Heisenberg & Euler (1936), Z. Phys. 98, 714
#   Raffelt & Stodolsky (1988), Phys. Rev. D 37, 1237
#   Dobrynina, Kartavtsev & Raffelt (2015), Phys. Rev. D 91, 083003
#   Meyer et al. (2021), arXiv:2108.02061 (gammaALPs implementation)

from utils.constants import G_TO_EV2

_ALPHA_FS      = 1.0 / 137.035999            # fine-structure constant
_B_CRIT_GAUSS  = 4.414e13                     # Schwinger critical field [G]
_B_CRIT_EV2    = _B_CRIT_GAUSS * G_TO_EV2     # critical field [eV^2]
_QED_COEFF     = _ALPHA_FS / (45.0 * np.pi)   # alpha/(45 pi)

# Birefringence coefficients for the two photon polarizations
_QED_PAR_FACTOR  = 3.5   # parallel photon  (7/2)
_QED_PERP_FACTOR = 2.0   # perpendicular photon  (4/2)


def _delta_qed_base(energy_ev, B_T_ev2):
    """
    Base QED vacuum birefringence term [eV] for transverse field magnitude B_T,
    before applying the parallel (3.5) or perpendicular (2.0) polarization
    factors.

    DQED = (alpha/45pi) * (B_T/B_crit)^2 * E * f(b),  b = B_T/B_crit
    f(b) = (1 + 1.2 b) / (1 + 1.33 b + 0.56 b^2)
    """
    if B_T_ev2 <= 0.0:
        return 0.0
    b    = B_T_ev2 / _B_CRIT_EV2
    corr = (1.0 + 1.2 * b) / (1.0 + 1.33 * b + 0.56 * b * b)
    return _QED_COEFF * (b ** 2) * energy_ev * corr


def _build_lab_matrix(energy_ev, Bx_ev2, By_ev2, g_ag_inv_ev, m_a_sq_ev2,
                      qed, extra_photon_diag):
    """
    Build the 3x3 mixing matrix in the lab frame (x, y, ALP) basis.

    The physics is defined in the B-aligned basis (par, perp, ALP):
        [ Dpar,  0,     Dag ]
        [ 0,     Dperp, 0   ]
        [ Dag,   0,     Da  ]
    where only the parallel photon couples to the ALP. This is rotated to the
    lab frame by the field angle psi = atan2(By, Bx).
    """
    B_T_ev2 = np.sqrt(Bx_ev2 * Bx_ev2 + By_ev2 * By_ev2)

    # Photon-ALP coupling uses the full transverse field magnitude
    Dag = 0.5 * g_ag_inv_ev * B_T_ev2
    Da  = m_a_sq_ev2 / (2.0 * energy_ev)

    if qed and B_T_ev2 > 0.0:
        DQED  = _delta_qed_base(energy_ev, B_T_ev2)
        Dpar  = _QED_PAR_FACTOR  * DQED
        Dperp = _QED_PERP_FACTOR * DQED
    else:
        Dpar  = 0.0
        Dperp = 0.0

    # B-aligned matrix: parallel photon along axis 0, perpendicular along axis 1
    M_aligned = np.array([
        [Dpar, 0.0,   Dag],
        [0.0,  Dperp, 0.0],
        [Dag,  0.0,   Da ]
    ], dtype=np.complex128)

    # Rotation angle of the transverse field in the lab frame
    if B_T_ev2 > 0.0:
        psi = np.arctan2(By_ev2, Bx_ev2)
    else:
        psi = 0.0
    c, s = np.cos(psi), np.sin(psi)
    R = np.array([
        [c, -s, 0.0],
        [s,  c, 0.0],
        [0.0, 0.0, 1.0]
    ], dtype=np.complex128)

    # Rotate to lab frame: M_lab = R · M_aligned · R^T
    M_lab = R @ M_aligned @ R.T

    # Add isotropic photon-diagonal term (absorption) after rotation
    if extra_photon_diag != 0.0:
        M_lab[0, 0] += extra_photon_diag
        M_lab[1, 1] += extra_photon_diag

    return M_lab


def get_mixing_matrix(energy_ev, Bx_ev2, By_ev2, g_ag_inv_ev, m_a_sq_ev2,
                      include_qed=True):
    """
    Constructs the 3x3 photon-ALP mixing matrix (M) in the lab frame.
    Units: [eV]

    Includes QED vacuum birefringence with the correct parallel (3.5) and
    perpendicular (2.0) polarization factors, implemented via the B-aligned
    basis and rotated to the lab frame.
    """
    return _build_lab_matrix(
        energy_ev, Bx_ev2, By_ev2, g_ag_inv_ev, m_a_sq_ev2,
        qed=include_qed, extra_photon_diag=0.0
    )


_unitarity_violation_count = 0


def get_mixing_matrix_with_absorption(energy_ev, Bx_ev2, By_ev2,
                                       g_ag_inv_ev, m_a_sq_ev2,
                                       gamma_absorb_ev, include_qed=True):
    """
    Mixing matrix including EBL absorption (non-Hermitian).
    """
    abs_diag = -0.5j * gamma_absorb_ev
    return _build_lab_matrix(
        energy_ev, Bx_ev2, By_ev2, g_ag_inv_ev, m_a_sq_ev2,
        qed=include_qed, extra_photon_diag=abs_diag
    )


def evolve_state(psi, M, dist_step_ev_inv):
    """
    Evolves the state vector psi via matrix exponentiation.
    M is Hermitian, so eigh is exact and returns real eigenvalues.
    """
    global _unitarity_violation_count
    try:
        norm_in = np.real(np.dot(psi.conj(), psi))
        w, v = np.linalg.eigh(M)
        phase = np.exp(-1j * w * dist_step_ev_inv)
        psi_out = v @ (phase * (v.conj().T @ psi))
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
