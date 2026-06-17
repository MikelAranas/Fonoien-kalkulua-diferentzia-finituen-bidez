import numpy as np
from numpy import linalg
from collections import defaultdict
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', '..', '..', 'modules'))
from svecs_module import get_svecs_multi, print_svecs_summary

M_C = 12.011
N   = 4

a1 = np.array([ 1.000000,  0.000000,  0.000000])
a2 = np.array([-0.500000,  np.sqrt(3)/2,  0.000000])
a3 = np.array([ 0.000000,  0.000000,  4.300000])

b1 = np.array([1.000000,  0.577350, 0.000000])
b2 = np.array([0.000000,  1.154701, 0.000000])
b3 = np.array([0.000000,  0.000000, 0.232558])

sc_positions = np.array([
    [0,      0,      0.0],
    [1/3,    2/3,    0.0],
    [1,      0,      0.0],
    [1+1/3,  2/3,    0.0],
    [2,      0,      0.0],
    [2+1/3,  2/3,    0.0],
    [3,      0,      0.0],
    [3+1/3,  2/3,    0.0],
    [0,      1,      0.0],
    [1/3,    1+2/3,  0.0],
    [1,      1,      0.0],
    [1+1/3,  1+2/3,  0.0],
    [2,      1,      0.0],
    [2+1/3,  1+2/3,  0.0],
    [3,      1,      0.0],
    [3+1/3,  1+2/3,  0.0],
    [0,      2,      0.0],
    [1/3,    2+2/3,  0.0],
    [1,      2,      0.0],
    [1+1/3,  2+2/3,  0.0],
    [2,      2,      0.0],
    [2+1/3,  2+2/3,  0.0],
    [3,      2,      0.0],
    [3+1/3,  2+2/3,  0.0],
    [0,      3,      0.0],
    [1/3,    3+2/3,  0.0],
    [1,      3,      0.0],
    [1+1/3,  3+2/3,  0.0],
    [2,      3,      0.0],
    [2+1/3,  3+2/3,  0.0],
    [3,      3,      0.0],
    [3+1/3,  3+2/3,  0.0],
], dtype=float)

n_atoms  = 32
n_basis  = 2
n_branch = 6

symmetry_points = [
    (r'$\Gamma$', np.array([0.000000,  0.000000,  0.0])),
    (r'$M$',      np.array([0.500000,  0.000000,  0.0])),
    (r'$K$',      np.array([0.333333,  0.333333,  0.0])),
    (r'$\Gamma$', np.array([0.000000,  0.000000,  0.0])),
]

q_path_file = 'path.dat'
IFC_file    = 'IFC_murriztue.dat'
output_plot = 'Grafenoa_DF_bandak.pdf'
output_dat  = 'Grafenoa_DF_bandak.dat'

factor = np.sqrt(13.605693 / 0.529177**2) * 15.633302 * 33.3564

A = np.vstack([a1, a2, a3])
B = np.vstack([b1, b2, b3])

IFC = np.loadtxt(IFC_file)
print("=" * 60)
print("DIAGNOSTICO IFC")
print("=" * 60)
print(f"  Shape IFC       : {IFC.shape}  (esperado: (6, 96))")
print(f"  IFC[0,0]  Phi_AA xx m=0 : {IFC[0,0]:.6f}")
print(f"  IFC[3,0]  Phi_BA xx m=0 : {IFC[3,0]:.6f}")
print(f"  Suma fila 0 (ASR, aprox 0): {np.sum(IFC[0,:]):.2e}")
print(f"  Suma fila 3 (ASR, aprox 0): {np.sum(IFC[3,:]):.2e}")

IFC_A = defaultdict(lambda: np.zeros((3, 3), dtype=float))
IFC_B = defaultdict(lambda: np.zeros((3, 3), dtype=float))

for m in range(n_atoms):
    for mu in range(3):
        for nu in range(3):
            IFC_A[m][mu, nu] = IFC[mu,     3*m + nu]
            IFC_B[m][mu, nu] = IFC[3 + mu, 3*m + nu]

basis_positions = np.array([
    [0.0,    0.0,    0.0],
    [1/3,    2/3,    0.0],
])
svecs, multi = get_svecs_multi(A, sc_positions, [N, N, 1], basis_positions=basis_positions)

print("\n" + "=" * 60)
print("DIAGNOSTICO SVECS")
print("=" * 60)
print_svecs_summary(svecs, multi)

q_path_red  = np.loadtxt(q_path_file)
QN          = q_path_red.shape[0]
q_path_cart = q_path_red @ B
s_path      = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_path_cart[i] - q_path_cart[i-1])

def phase_avg(q_red, svecs, multi, k, s):
    m_val, adrs = multi[k, s]
    return sum(
        np.exp(1j * 2 * np.pi * np.dot(q_red, svecs[adrs + ll]))
        for ll in range(m_val)
    ) / m_val

def dynamical_matrix(q):
    q_red = q
    D = np.zeros((6, 6), dtype=complex)

    for m in range(n_atoms):
        ph_A = phase_avg(q_red, svecs, multi, m, s=0)
        ph_B = phase_avg(q_red, svecs, multi, m, s=1)

        if m % 2 == 0:
            D[0:3, 0:3] += IFC_A[m] * ph_A
            D[3:6, 0:3] += IFC_B[m] * ph_B
        else:
            D[0:3, 3:6] += IFC_A[m] * ph_A
            D[3:6, 3:6] += IFC_B[m] * ph_B

    D = (D + D.conj().T) / 2.0
    return D / M_C

def phonon_frequencies(D):
    w2, _ = np.linalg.eigh(D)
    w2[w2 < 0] = 0.0
    return np.sqrt(w2) * factor

print("\n" + "=" * 60)
print("DIAGNOSTICO EN GAMMA")
print("=" * 60)
D_gamma = dynamical_matrix(np.array([0.0, 0.0, 0.0]))
freqs_gamma = phonon_frequencies(D_gamma)
print(f"  Frecuencias en Gamma: {np.round(freqs_gamma, 4)} cm-1")
print(f"  (Esperado: 3 acusticos ~0, ZO ~873, LO/TO ~1559 cm-1)")

print("\n" + "=" * 60)
print("DIAGNOSTICO EN K")
print("=" * 60)
D_K = dynamical_matrix(np.array([1/3, 1/3, 0.0]))
freqs_K = phonon_frequencies(D_K)
print(f"  Frecuencias en K: {np.round(freqs_K, 4)} cm-1")
print(f"  (Esperado: pares degenerados [0,0], [1093,1093], [1329,1329])")

omega = np.zeros((QN, n_branch))
for i in range(QN):
    omega[i] = phonon_frequencies(dynamical_matrix(q_path_red[i]))

bandas = np.column_stack((s_path, omega))

def s_at_q_on_polyline(q_path, s_path, q_target):
    best_s, best_d2 = None, np.inf
    for i in range(len(q_path) - 1):
        a, b = q_path[i], q_path[i+1]
        v    = b - a
        vv   = np.dot(v, v)
        if vv == 0.0:
            continue
        t    = np.clip(np.dot(q_target - a, v) / vv, 0.0, 1.0)
        d2   = np.dot(q_target - (a + t*v), q_target - (a + t*v))
        if d2 < best_d2:
            best_d2 = d2
            best_s  = s_path[i] + t * linalg.norm(v)
    return best_s

tick_labels, tick_positions = [], []
for label, qpt in symmetry_points:
    tick_labels.append(label)
    tick_positions.append(s_at_q_on_polyline(q_path_cart, s_path, qpt @ B))

plt.figure(figsize=(6, 7))
for b in range(n_branch):
    plt.plot(s_path, omega[:, b], color='royalblue', lw=1.5)
for s_tick in tick_positions:
    plt.axvline(s_tick, color='0.5', ls='--', lw=0.8)
plt.axhline(0, color='k', lw=0.5)
plt.xlim(s_path[0], s_path[-1])
plt.ylim(bottom=0)
plt.ylabel(r'$\omega$ (cm$^{-1}$)')
plt.xticks(tick_positions, tick_labels, fontsize=13)
plt.title('Grafenoa DF 4x4x1')
plt.tight_layout()
plt.savefig(output_plot)
print(f"\n  Grafica guardada en {output_plot}")
np.savetxt(output_dat, bandas)
