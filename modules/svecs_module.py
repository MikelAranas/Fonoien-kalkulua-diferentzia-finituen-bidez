"""
svecs_module.py
===============
Shortest periodic images for pair displacements in a supercell.

Material- and shape-agnostic:
  * monoatomic (Al-style) and polyatomic (Si, Grafenoa) — via `basis_positions`
  * diagonal N×N×N, anisotropic N1×N2×N3 (e.g., graphene 2D), and arbitrary
    non-diagonal supercells — via `supercell_matrix`

Public API
----------
    get_svecs_multi(prim_vecs, sc_positions, supercell_matrix,
                    basis_positions=None)
        → (svecs, multi)

    print_svecs_summary(svecs, multi)

Usage from every material's D.py
--------------------------------
    import sys, os
    sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                     '..', '..', '..', '..', 'modules'))
    from svecs_module import get_svecs_multi

    # Al (monoatomic FCC, 6×6×6):
    svecs, multi = get_svecs_multi(prim, sc_positions, 6)

    # Si (diamond, 2-atom basis, 6×6×6):
    basis = np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]])
    svecs, multi = get_svecs_multi(prim, sc_positions, 6,
                                   basis_positions=basis)

    # Grafenoa (honeycomb 2D, 7×7×1):
    basis = np.array([[0.0, 0.0, 0.0], [1/3, 2/3, 0.0]])
    svecs, multi = get_svecs_multi(prim, sc_positions, [7, 7, 1],
                                   basis_positions=basis)

    # Future non-diagonal supercell:
    S = np.array([[2, 1, 0], [0, 1, 0], [0, 0, 1]], dtype=int)
    svecs, multi = get_svecs_multi(prim, sc_positions, S,
                                   basis_positions=basis)

Indexing convention
-------------------
For each supercell atom j ∈ [0, N_sc) and each basis atom s ∈ [0, n_basis):
    m_val = multi[j, s, 0]   # number of equidistant shortest images
    adrs  = multi[j, s, 1]   # starting index in svecs
    images = svecs[adrs : adrs + m_val]
Each image is a 3-vector in PRIMITIVE FRACTIONAL coordinates, pointing from
basis atom s to supercell atom j (modulo supercell periodicity).

Atom ordering: sc_positions MUST be in the same order as the ATOMIC_POSITIONS
block of the QE input that generated the IFC matrix. Otherwise phases are
scrambled.
"""

import numpy as np
from itertools import product


# ─────────────────────────────────────────────────────────────────────────────
#  Private helpers
# ─────────────────────────────────────────────────────────────────────────────

def _metric_tensor(prim_vecs):
    """G_ij = a_i · a_j. Lets us compute |R|² in primitive fractional coords."""
    A = np.asarray(prim_vecs, dtype=float)
    return A @ A.T


def _to_supercell_matrix(supercell_matrix):
    """
    Normalise the supercell_matrix argument to a (3, 3) integer numpy array.

    Accepts:
        scalar N           → N · I  (diagonal N×N×N)
        [N1, N2, N3]       → diag(N1, N2, N3)  (anisotropic diagonal)
        (3, 3) int matrix  → returned as-is (may be non-diagonal)
    """
    S = np.asarray(supercell_matrix)
    if S.ndim == 0:
        return int(S) * np.eye(3, dtype=int)
    if S.ndim == 1 and len(S) == 3:
        return np.diag(np.asarray(S, dtype=int))
    if S.shape == (3, 3):
        return S.astype(int)
    raise ValueError(
        f"supercell_matrix must be scalar, length-3 sequence, or 3×3 matrix; "
        f"got array of shape {S.shape}"
    )


def _shortest_images(diff, S, G, search_range, tol):
    """
    Find all periodic images of `diff` with minimum |R|² under the
    supercell periodicity defined by S.

    Periodic shifts: R = diff + n1·S[0] + n2·S[1] + n3·S[2]
    for n1, n2, n3 ∈ [-search_range, +search_range].

    `tol` is the distance² tolerance for declaring two images "equidistant"
    (→ multiplicity > 1; this happens at the Wigner–Seitz boundary).
    """
    candidates = []
    for n1, n2, n3 in product(range(-search_range, search_range + 1), repeat=3):
        R = diff + n1 * S[0] + n2 * S[1] + n3 * S[2]
        dist2 = float(R @ G @ R)
        candidates.append((dist2, R))
    min_d2 = min(c[0] for c in candidates)
    return [c[1] for c in candidates if abs(c[0] - min_d2) < tol]


# ─────────────────────────────────────────────────────────────────────────────
#  Public API
# ─────────────────────────────────────────────────────────────────────────────

def get_svecs_multi(prim_vecs, sc_positions, supercell_matrix,
                    basis_positions=None, search_range=3, tol=1e-6):
    """
    Compute shortest periodic images for every (supercell atom, basis atom) pair.

    Parameters
    ----------
    prim_vecs : (3, 3) array_like
        Primitive lattice vectors as rows. Geometric scale (alat / Å / Bohr)
        is irrelevant — only relative geometry matters.
    sc_positions : (N_sc, 3) array_like
        Supercell atom positions in PRIMITIVE fractional coordinates, ordered
        to match the ATOMIC_POSITIONS block of the QE input.
    supercell_matrix : scalar | length-3 sequence | (3, 3) int matrix
        Defines supercell periodicity. See `_to_supercell_matrix` for accepted
        forms; supports diagonal and non-diagonal cases uniformly.
    basis_positions : (n_basis, 3) array_like, optional
        Primitive-cell basis-atom positions in fractional coordinates.
        Default [[0, 0, 0]] for monoatomic.
    search_range : int
        Number of supercell cells to enumerate in each direction when seeking
        shortest images. 3 is sufficient for typical diagonal supercells up
        to ~7. Increase for larger or strongly anisotropic supercells.
    tol : float
        Distance² tolerance for equidistant-image detection.

    Returns
    -------
    svecs : (total_vecs, 3) ndarray
        All shortest-image vectors concatenated. Primitive fractional coords.
    multi : (N_sc, n_basis, 2) int ndarray
        multi[j, s, 0] = multiplicity of the (sc atom j, basis atom s) pair
        multi[j, s, 1] = starting index in `svecs` for that pair
    """
    if basis_positions is None:
        basis_positions = np.array([[0.0, 0.0, 0.0]])
    basis_positions = np.asarray(basis_positions, dtype=float)
    sc_positions    = np.asarray(sc_positions, dtype=float)

    S = _to_supercell_matrix(supercell_matrix)
    G = _metric_tensor(prim_vecs)

    N_sc    = len(sc_positions)
    n_basis = len(basis_positions)

    svecs_buf = []
    multi = np.zeros((N_sc, n_basis, 2), dtype=int)
    idx = 0
    for j in range(N_sc):
        for s in range(n_basis):
            diff   = sc_positions[j] - basis_positions[s]
            images = _shortest_images(diff, S, G, search_range, tol)
            multi[j, s, 0] = len(images)
            multi[j, s, 1] = idx
            svecs_buf.extend(images)
            idx += len(images)

    return np.array(svecs_buf), multi


def print_svecs_summary(svecs, multi):
    """
    Human-readable dump of svecs and multiplicities. Useful for verifying that
    sc_positions are in the same order as the QE input.
    """
    N_sc, n_basis, _ = multi.shape
    print(f"{'SC':>4}  {'s':>2}  {'mult':>4}  Shortest images "
          f"(primitive frac coords)")
    print("-" * 70)
    for j in range(N_sc):
        for s in range(n_basis):
            m, adrs = multi[j, s]
            vecs    = svecs[adrs:adrs + m]
            txt     = ",  ".join(
                f"[{v[0]:+.3f},{v[1]:+.3f},{v[2]:+.3f}]" for v in vecs
            )
            print(f"  {j:>3}  {s:>2}  {m:>4}  {txt}")


# ─────────────────────────────────────────────────────────────────────────────
#  Self-test (only runs when invoked as a script)
# ─────────────────────────────────────────────────────────────────────────────

if __name__ == "__main__":
    # Al FCC 2×2×2, monoatomic. sc_positions in the QE convention
    # (ix slowest, iz fastest).
    prim_Al = np.array([
        [-0.707107,  0.000000,  0.707107],
        [ 0.000000,  0.707107,  0.707107],
        [-0.707107,  0.707107,  0.000000],
    ])
    sc_Al = np.array(
        [[ix, iy, iz] for ix in range(2) for iy in range(2) for iz in range(2)],
        dtype=float,
    )
    svecs, multi = get_svecs_multi(prim_Al, sc_Al, 2)
    print("Al FCC 2×2×2 (monoatomic):")
    print_svecs_summary(svecs, multi)

    # Graphene 2×2×1, honeycomb 2-atom basis. supercell_matrix = [2, 2, 1].
    print()
    print("Graphene 2×2×1 (honeycomb):")
    prim_C = np.array([
        [ 1.0,            0.0,           0.0],
        [-0.5,            np.sqrt(3)/2,  0.0],
        [ 0.0,            0.0,           4.3],
    ])
    basis_C = np.array([[0.0, 0.0, 0.0], [1/3, 2/3, 0.0]])
    sc_C = []
    for n2 in range(2):
        for n1 in range(2):
            sc_C.append([n1,        n2,        0.0])
            sc_C.append([n1 + 1/3,  n2 + 2/3,  0.0])
    sc_C = np.array(sc_C, dtype=float)
    svecs, multi = get_svecs_multi(prim_C, sc_C, [2, 2, 1],
                                   basis_positions=basis_C)
    print_svecs_summary(svecs, multi)
