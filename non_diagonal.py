"""
Non-diagonal supercell generator for finite-difference phonon calculations.

Implements the algorithm from:
    Lloyd-Williams & Monserrat, Phys. Rev. B 92, 184301 (2015)
    "Lattice dynamics and electron-phonon coupling calculations
     using non-diagonal supercells"

Given a q-point with reduced fractional coordinates (m1/n1, m2/n2, m3/n3),
this module returns an integer 3x3 supercell matrix S such that:

    A_super = S @ A_prim     (rows are lattice vectors)

contains LCM(n1, n2, n3) primitive cells and is commensurate with q.
The supercell shape is then optimized (made compact) by unimodular row
operations, which do not change the superlattice — only its basis vectors.

Usage
-----
>>> import numpy as np
>>> from nondiagonal_supercell import build_supercell
>>>
>>> # Primitive lattice (rows = lattice vectors). Example: simple cubic, a=1.
>>> A_prim = np.eye(3)
>>>
>>> # q-point as fractions of the reciprocal primitive basis.
>>> # Format: list of (numerator, denominator) pairs, NOT already reduced.
>>> q = [(1, 4), (1, 2), (1, 4)]
>>>
>>> S = build_supercell(q)
>>> print("Supercell matrix S =\\n", S)
>>> print("Number of primitive cells in supercell:", int(round(np.linalg.det(S))))
>>>
>>> A_super = S @ A_prim
>>> print("Supercell lattice vectors:\\n", A_super)
"""

from __future__ import annotations
from math import gcd
from itertools import product
import numpy as np


# ------------------------------------------------------------------
# Number-theory helpers
# ------------------------------------------------------------------

def _lcm(a: int, b: int) -> int:
    """Least common multiple of two positive integers."""
    return a * b // gcd(a, b)


def lcm3(n1: int, n2: int, n3: int) -> int:
    """LCM of three positive integers."""
    return _lcm(_lcm(n1, n2), n3)


def reduce_fraction(m: int, n: int) -> tuple[int, int]:
    """Reduce m/n to lowest terms, with n > 0.

    Handles negative numerators correctly. The wave vector is defined modulo
    integers, so we also wrap m into [0, n) for canonical form.
    """
    if n <= 0:
        raise ValueError(f"Denominator must be positive, got n = {n}")
    m = m % n               # wrap into [0, n)
    if m == 0:
        return 0, 1
    g = gcd(m, n)
    return m // g, n // g


# ------------------------------------------------------------------
# Hermite Normal Form enumeration (eqs. 7-12 of the paper)
# ------------------------------------------------------------------

def enumerate_hnf_supercells(q_frac):
    """Yield all valid HNF supercell matrices commensurate with q.

    Parameters
    ----------
    q_frac : sequence of (m_i, n_i) pairs
        The q-point as fractions of the reciprocal primitive basis.
        Internally reduced to lowest terms first.

    Yields
    ------
    S : (3, 3) integer ndarray
        Upper-triangular HNF supercell matrix with det(S) = LCM(n1, n2, n3).

    Notes
    -----
    The HNF form is:

        [ S11  S12  S13 ]
        [  0   S22  S23 ]
        [  0    0   S33 ]

    with 0 <= S12 < S22 and 0 <= S13, S23 < S33.

    Commensurability with q means S @ q_p must have all-integer components
    (where q_p is q in reciprocal primitive coordinates). For a wave vector
    (m1/n1, m2/n2, m3/n3) this is equivalent to:

        S11*m1/n1 + S12*m2/n2 + S13*m3/n3   in Z   (row 1)
        S22*m2/n2 + S23*m3/n3               in Z   (row 2)
        S33*m3/n3                           in Z   (row 3)

    We enumerate all (S11, S22, S33) triples with S11*S22*S33 = LCM(n1,n2,n3)
    and then test each off-diagonal triple (S12, S13, S23).
    """
    # Reduce each m_i/n_i to lowest terms.
    reduced = [reduce_fraction(m, n) for m, n in q_frac]
    m1, n1 = reduced[0]
    m2, n2 = reduced[1]
    m3, n3 = reduced[2]

    L = lcm3(n1, n2, n3)

    # Loop over all integer factorizations S11*S22*S33 = L.
    # For 3D this is a tiny search space; brute force is fine.
    for S11 in range(1, L + 1):
        if L % S11 != 0:
            continue
        rem1 = L // S11
        for S22 in range(1, rem1 + 1):
            if rem1 % S22 != 0:
                continue
            S33 = rem1 // S22

            # Test all valid off-diagonal entries.
            for S12, S13, S23 in product(range(S22), range(S33), range(S33)):
                # Row 3: S33*m3 must be divisible by n3.
                if (S33 * m3) % n3 != 0:
                    continue

                # Row 2: S22*m2/n2 + S23*m3/n3 must be an integer.
                # Combine over common denominator lcm(n2, n3).
                num = S22 * m2 * (lcm3(n2, n3, 1) // n2) \
                      + S23 * m3 * (lcm3(n2, n3, 1) // n3)
                if num % lcm3(n2, n3, 1) != 0:
                    continue

                # Row 1: similarly over lcm(n1, n2, n3) = L.
                num1 = S11 * m1 * (L // n1) \
                       + S12 * m2 * (L // n2) \
                       + S13 * m3 * (L // n3)
                if num1 % L != 0:
                    continue

                S = np.array([[S11, S12, S13],
                              [0,   S22, S23],
                              [0,   0,   S33]], dtype=int)
                yield S


# ------------------------------------------------------------------
# Compactness scoring & lattice reduction
# ------------------------------------------------------------------

def compactness_score(A: np.ndarray) -> float:
    """Lower is more compact.

    Heuristic: sum of squared lengths of the three lattice vectors (rows of A).
    For a fixed volume |det A|, this is minimized when the cell is cube-like.
    Equivalent to the Frobenius norm squared.
    """
    return float(np.sum(A * A))


def reduce_supercell(S: np.ndarray, A_prim: np.ndarray,
                     max_passes: int = 50) -> np.ndarray:
    """Apply unimodular row operations to make the supercell shape compact.

    We iteratively try replacing each row by (row_i + k * row_j) for small
    integer k, accepting the change if it shortens the resulting lattice
    vector. This is a Minkowski-style 3D reduction — simpler than LLL,
    plenty good enough for crystal supercells.

    The superlattice itself is invariant under these operations; only the
    basis vectors change.

    Parameters
    ----------
    S : (3, 3) int ndarray
        Supercell matrix (rows are integer combinations of primitive vectors).
    A_prim : (3, 3) float ndarray
        Primitive lattice (rows are Cartesian lattice vectors).
    max_passes : int
        Safety bound on the reduction loop.

    Returns
    -------
    S_reduced : (3, 3) int ndarray
        An equivalent supercell matrix with shorter lattice vectors.
    """
    S = S.copy()
    for _ in range(max_passes):
        A_super = S @ A_prim
        improved = False
        for i in range(3):
            for j in range(3):
                if i == j:
                    continue
                # Find the best integer k for row_i <- row_i + k * row_j.
                vi = A_super[i]
                vj = A_super[j]
                vj2 = float(vj @ vj)
                if vj2 == 0:
                    continue
                # Closed-form optimum: k* = -(vi . vj) / (vj . vj).
                k_opt = -int(round((vi @ vj) / vj2))
                if k_opt == 0:
                    continue
                new_vi = vi + k_opt * vj
                if float(new_vi @ new_vi) < float(vi @ vi) - 1e-12:
                    S[i] = S[i] + k_opt * S[j]
                    A_super = S @ A_prim
                    improved = True
        if not improved:
            break

    # Ensure positive orientation (right-handed): flip the last row if needed.
    if np.linalg.det(S @ A_prim) < 0:
        S[2] = -S[2]
    return S


# ------------------------------------------------------------------
# Main entry point
# ------------------------------------------------------------------

def build_supercell(q_frac, A_prim=None, return_all=False):
    """Find the optimal non-diagonal supercell for a given q-point.

    Parameters
    ----------
    q_frac : sequence of (m_i, n_i) pairs, length 3
        The q-point in fractional coordinates of the reciprocal primitive
        basis. E.g. q = [(1, 4), (1, 2), (1, 4)] for (1/4, 1/2, 1/4).
    A_prim : (3, 3) ndarray or None
        Cartesian primitive lattice (rows = lattice vectors). Used to score
        compactness. If None, the identity is used (cubic with a=1), which
        is fine if you only care about |S|.
    return_all : bool
        If True, also return the list of all candidate HNF matrices considered.

    Returns
    -------
    S_best : (3, 3) int ndarray
        The supercell matrix with the most compact shape.
        det(S_best) = LCM(n1, n2, n3) primitive cells.
    candidates : list of (3, 3) int ndarrays   (only if return_all=True)
    """
    if A_prim is None:
        A_prim = np.eye(3)
    A_prim = np.asarray(A_prim, dtype=float)

    candidates = []
    best_S = None
    best_score = np.inf

    for S_hnf in enumerate_hnf_supercells(q_frac):
        S_red = reduce_supercell(S_hnf, A_prim)
        score = compactness_score(S_red @ A_prim)
        candidates.append(S_red)
        if score < best_score:
            best_score = score
            best_S = S_red

    if best_S is None:
        raise RuntimeError(f"No commensurate supercell found for q = {q_frac}")

    if return_all:
        return best_S, candidates
    return best_S


# ------------------------------------------------------------------
# Verification helper
# ------------------------------------------------------------------

def check_commensurate(S: np.ndarray, q_frac) -> bool:
    """Verify that the wave vector q is commensurate with supercell S.

    Commensurate means S @ q has all-integer components.
    """
    q = np.array([m / n for m, n in q_frac], dtype=float)
    S_q = S @ q
    return np.all(np.abs(S_q - np.round(S_q)) < 1e-10)


if __name__ == "__main__":
    # ----------------------------------------------------------------
    # Demonstration: the 2D example from Fig. 1 of the paper, embedded
    # in 3D by adding a trivial third dimension. Square lattice, q = (1/2, 1/2).
    # ----------------------------------------------------------------
    print("=" * 70)
    print("Example 1: q = (1/2, 1/2, 0) on a simple cubic lattice")
    print("           (the 2D Fig. 1 example, with a trivial 3rd axis)")
    print("=" * 70)
    A_prim = np.eye(3)
    q = [(1, 2), (1, 2), (0, 1)]
    S = build_supercell(q, A_prim)
    print(f"  Reduced q       = {[reduce_fraction(m, n) for m, n in q]}")
    print(f"  LCM(2, 2, 1)    = {lcm3(2, 2, 1)}  -> supercell has 2 prim cells")
    print(f"  Supercell S     =\n{S}")
    print(f"  |det S|         = {int(round(abs(np.linalg.det(S))))}")
    print(f"  Commensurate?   = {check_commensurate(S, q)}")

    print()
    print("=" * 70)
    print("Example 2: q = (1/4, 1/2, 1/4) — the diamond worst case from")
    print("           Section IV.B of the paper.")
    print("=" * 70)
    q = [(1, 4), (1, 2), (1, 4)]
    S = build_supercell(q, A_prim)
    print(f"  Reduced q       = {[reduce_fraction(m, n) for m, n in q]}")
    print(f"  LCM(4, 2, 4)    = {lcm3(4, 2, 4)}  -> supercell has 4 prim cells")
    print(f"  Diagonal would need {4*2*4} = 32 primitive cells.")
    print(f"  Supercell S     =\n{S}")
    print(f"  |det S|         = {int(round(abs(np.linalg.det(S))))}")
    print(f"  Commensurate?   = {check_commensurate(S, q)}")

    print()
    print("=" * 70)
    print("Example 3: scan a 4x4x4 q-grid (no symmetry) — show sizes needed")
    print("=" * 70)
    N = 4
    sizes_diag = []
    sizes_ndsc = []
    for i, j, k in product(range(N), range(N), range(N)):
        if i == j == k == 0:
            continue  # Gamma point, trivially the primitive cell
        q = [(i, N), (j, N), (k, N)]
        reduced = [reduce_fraction(m, n) for m, n in q]
        ns = [n for _, n in reduced]
        size_diag = ns[0] * ns[1] * ns[2]
        size_ndsc = lcm3(*ns)
        sizes_diag.append(size_diag)
        sizes_ndsc.append(size_ndsc)
    print(f"  Across {len(sizes_diag)} non-trivial q-points on a {N}x{N}x{N} grid:")
    print(f"    Largest diagonal supercell needed:     {max(sizes_diag)} prim cells")
    print(f"    Largest non-diagonal supercell needed: {max(sizes_ndsc)} prim cells")
    print(f"    Sum of cell sizes (rough cost proxy):")
    print(f"      Diagonal:      {sum(sizes_diag)}")
    print(f"      Non-diagonal:  {sum(sizes_ndsc)}")
    print(f"      Ratio:         {sum(sizes_diag) / sum(sizes_ndsc):.2f}x")
