import numpy as np
from numpy import linalg
from collections import defaultdict
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', '..', '..', 'modules'))
from svecs_module import get_svecs_multi, print_svecs_summary

M_Si = 28.086
N    = 2

a1 = np.array([-0.707107,  0.000000,  0.707107])
a2 = np.array([ 0.000000,  0.707107,  0.707107])
a3 = np.array([-0.707107,  0.707107,  0.000000])

b1 = np.array([-0.707107, -0.707107,  0.707107])
b2 = np.array([ 0.707107,  0.707107,  0.707107])
b3 = np.array([-0.707107,  0.707107, -0.707107])

sc_positions = np.array([
    [0,      0,      0    ],
    [0.25,   0.25,   0.25 ],
    [1,      0,      0    ],
    [1.25,   0.25,   0.25 ],
    [0,      1,      0    ],
    [0.25,   1.25,   0.25 ],
    [0,      0,      1    ],
    [0.25,   0.25,   1.25 ],
    [1,      1,      0    ],
    [1.25,   1.25,   0.25 ],
    [1,      0,      1    ],
    [1.25,   0.25,   1.25 ],
    [0,      1,      1    ],
    [0.25,   1.25,   1.25 ],
    [1,      1,      1    ],
    [1.25,   1.25,   1.25 ],
], dtype=float)

n_atoms  = N**3 * 2
n_basis  = 2
n_branch = 3 * n_basis

symmetry_points = [
    (r'$W$',   np.array([0.75,   0.50,  0.25 ])),
    (r'$X$',   np.array([0.50,   0.50,  0.00 ])),
    (r'$\Gamma$', np.array([0.00,   0.00,  0.00 ])),
    (r'$L$',   np.array([0.50,   0.50,  0.50 ])),
    (r'$K$',   np.array([0.75,   0.375, 0.375])),
]

q_path_file = 'path.dat'
IFC_file    = 'IFC_murriztue.dat'
output_plot = 'Si_DF_bandak.pdf'
output_dat  = 'Si_DF_bandak.dat'

factor = np.sqrt(13.605693 / 0.529177**2) * 15.633302

A = np.vstack([a1, a2, a3])
B = np.vstack([b1, b2, b3])

IFC = np.loadtxt(IFC_file)
print("=" * 60)
print("DIAGNOSTICO IFC")
print("=" * 60)
print(f"  Shape IFC       : {IFC.shape}  (esperado: (6, 48))")
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

print("\n  IFC_A[0]  (A desplazado → atomo 0, subred A):")
print("  " + str(np.round(IFC_A[0], 5)).replace("\n", "\n  "))
print("\n  IFC_A[1]  (A desplazado → atomo 1, subred B):")
print("  " + str(np.round(IFC_A[1], 5)).replace("\n", "\n  "))
print("\n  IFC_B[0]  (B desplazado → atomo 0, subred A):")
print("  " + str(np.round(IFC_B[0], 5)).replace("\n", "\n  "))
print("\n  IFC_B[1]  (B desplazado → atomo 1, subred B):")
print("  " + str(np.round(IFC_B[1], 5)).replace("\n", "\n  "))

basis_positions = np.array([
    [0.00, 0.00, 0.00],
    [0.25, 0.25, 0.25],
])
svecs, multi = get_svecs_multi(A, sc_positions, N, basis_positions=basis_positions)

print("\n" + "=" * 60)
print("DIAGNOSTICO SVECS")
print("=" * 60)
print_svecs_summary(svecs, multi)

q_path_red  = np.loadtxt(q_path_file, unpack=False, skiprows=0, usecols=(0, 1, 2))
QN          = q_path_red.shape[0]
q_path_cart = q_path_red @ B
s_path      = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_path_cart[i] - q_path_cart[i-1])
print(f"\n  Q-path: {QN} puntos")
print(f"  Primer q (cristalino): {q_path_red[0]}")
print(f"  Último  q (cristalino): {q_path_red[-1]}")

def phase_avg(q_red, svecs, multi, k, s):
    m_val, adrs = multi[k, s]
    return sum(
        np.exp(1j * 2 * np.pi * np.dot(q_red, svecs[adrs + ll]))
        for ll in range(m_val)
    ) / m_val

def dynamical_matrix(q, B):
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
    return D / M_Si

def phonon_frequencies(D):
    w2, _ = np.linalg.eigh(D)
    w2[w2 < 0] = 0.0
    return np.sqrt(w2) * factor

print("\n" + "=" * 60)
print("CALCULANDO BANDAS...")
print("=" * 60)

omega = np.zeros((QN, n_branch))
for i in range(QN):
    omega[i] = phonon_frequencies(dynamical_matrix(q_path_red[i], B))

n_neg = np.sum(np.any(omega == 0.0, axis=1))
print(f"  Puntos q con algun modo forzado a 0: {n_neg}")
print(f"  Frecuencia maxima : {omega.max():.4f} THz")
print(f"  Frecuencia minima no nula: {omega[omega > 1.0].min():.4f} THz")

bandas = np.column_stack((s_path, omega))

def s_at_q_on_polyline(q_path, s_path, q_target):
    best_s, best_d2 = None, np.inf
    for i in range(len(q_path) - 1):
        a, b = q_path[i], q_path[i + 1]
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

plt.figure(figsize=(7.5, 4.5))
for b in range(n_branch):
    plt.plot(s_path, omega[:, b], color='C0', lw=1.5)
for s_tick in tick_positions:
    plt.axvline(s_tick, color='0.7', ls=':', lw=0.8, zorder=0)

plt.xlim(s_path[0], s_path[-1])
plt.ylim(bottom=0)
plt.ylabel(r'$\omega$ (THz)')
plt.xticks(tick_positions, tick_labels)
plt.tight_layout()
plt.savefig(output_plot)
print(f"\n  Grafica guardada en {output_plot}")

np.savetxt(output_dat, bandas)
print(f"  Datos guardados en {output_dat}")
