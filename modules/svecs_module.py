import numpy as np
from itertools import product

def _metric_tensor(prim_vecs):
    A = np.asarray(prim_vecs, dtype=float)
    return A @ A.T

def _to_supercell_matrix(supercell_matrix):
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
    candidates = []
    for n1, n2, n3 in product(range(-search_range, search_range + 1), repeat=3):
        R = diff + n1 * S[0] + n2 * S[1] + n3 * S[2]
        dist2 = float(R @ G @ R)
        candidates.append((dist2, R))
    min_d2 = min(c[0] for c in candidates)
    return [c[1] for c in candidates if abs(c[0] - min_d2) < tol]

def get_svecs_multi(prim_vecs, sc_positions, supercell_matrix,
                    basis_positions=None, search_range=3, tol=1e-6):
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

if __name__ == "__main__":
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
