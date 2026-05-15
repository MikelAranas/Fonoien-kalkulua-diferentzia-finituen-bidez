import numpy as np
from numpy import linalg
from collections import defaultdict
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', '..', '..', 'modules'))
from svecs_module import get_svecs_multi, print_svecs_summary

# ═════════════════════════════════════════════════════════════════
#  PARAMETROS DEL MATERIAL
# ═════════════════════════════════════════════════════════════════

M_Si = 28.086   # masa atómica del Si (uma)
N    = 2        # supercelda 2×2×2

a1 = np.array([-0.707107,  0.000000,  0.707107])
a2 = np.array([ 0.000000,  0.707107,  0.707107])
a3 = np.array([-0.707107,  0.707107,  0.000000])

b1 = np.array([-0.707107, -0.707107,  0.707107])
b2 = np.array([ 0.707107,  0.707107,  0.707107])
b3 = np.array([-0.707107,  0.707107, -0.707107])
#b1 = np.array([-1.0, -1.0,  1.0])
#b2 = np.array([ 1.0,  1.0,  1.0])
#b3 = np.array([-1.0,  1.0, -1.0])

sc_positions = np.array([
    [0,      0,      0    ],  # 0  subred A
    [0.25,   0.25,   0.25 ],  # 1  subred B
    [1,      0,      0    ],  # 2  A
    [1.25,   0.25,   0.25 ],  # 3  B
    [0,      1,      0    ],  # 4  A
    [0.25,   1.25,   0.25 ],  # 5  B
    [0,      0,      1    ],  # 6  A
    [0.25,   0.25,   1.25 ],  # 7  B
    [1,      1,      0    ],  # 8  A
    [1.25,   1.25,   0.25 ],  # 9  B
    [1,      0,      1    ],  # 10 A
    [1.25,   0.25,   1.25 ],  # 11 B
    [0,      1,      1    ],  # 12 A
    [0.25,   1.25,   1.25 ],  # 13 B
    [1,      1,      1    ],  # 14 A
    [1.25,   1.25,   1.25 ],  # 15 B
], dtype=float)

n_atoms  = N**3 * 2   # 16
n_basis  = 2
n_branch = 3 * n_basis  # 6 ramas

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

factor = np.sqrt(13.605693 / 0.529177**2) * 15.633302* 33.3564

# ═════════════════════════════════════════════════════════════════
#  CARGA Y DIAGNOSTICO DE IFC
# ═════════════════════════════════════════════════════════════════

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

# Bloques 3×3 por atomo y por atomo desplazado
IFC_A = defaultdict(lambda: np.zeros((3, 3), dtype=float))
IFC_B = defaultdict(lambda: np.zeros((3, 3), dtype=float))

for m in range(n_atoms):
    for mu in range(3):
        for nu in range(3):
            IFC_A[m][mu, nu] = IFC[mu,     3*m + nu]  # atomo A desplazado
            IFC_B[m][mu, nu] = IFC[3 + mu, 3*m + nu]  # atomo B desplazado

print("\n  IFC_A[0]  (A desplazado → atomo 0, subred A):")
print("  " + str(np.round(IFC_A[0], 5)).replace("\n", "\n  "))
print("\n  IFC_A[1]  (A desplazado → atomo 1, subred B):")
print("  " + str(np.round(IFC_A[1], 5)).replace("\n", "\n  "))
print("\n  IFC_B[0]  (B desplazado → atomo 0, subred A):")
print("  " + str(np.round(IFC_B[0], 5)).replace("\n", "\n  "))
print("\n  IFC_B[1]  (B desplazado → atomo 1, subred B):")
print("  " + str(np.round(IFC_B[1], 5)).replace("\n", "\n  "))

# ═════════════════════════════════════════════════════════════════
#  SVECS Y MULTI
# ═════════════════════════════════════════════════════════════════

basis_positions = np.array([
    [0.00, 0.00, 0.00],   # Si_A
    [0.25, 0.25, 0.25],   # Si_B
])
svecs, multi = get_svecs_multi(A, sc_positions, N, basis_positions=basis_positions)

print("\n" + "=" * 60)
print("DIAGNOSTICO SVECS")
print("=" * 60)
print_svecs_summary(svecs, multi)

# ═════════════════════════════════════════════════════════════════
#  Q-PATH
# ═════════════════════════════════════════════════════════════════

q_path_red  = np.loadtxt(q_path_file, unpack=False, skiprows=0, usecols=(0, 1, 2))
QN          = q_path_red.shape[0]
q_path_cart = q_path_red @ B
s_path      = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_path_cart[i] - q_path_cart[i-1])
print(f"\n  Q-path: {QN} puntos")
print(f"  Primer q (cristalino): {q_path_red[0]}")
print(f"  Último  q (cristalino): {q_path_red[-1]}")

# ═════════════════════════════════════════════════════════════════
#  MATRIZ DINAMICA 6×6
# ═════════════════════════════════════════════════════════════════

def phase_avg(q_red, svecs, multi, k, s):
    """Promedio de fase sobre vectores equivalentes del par (atomo k, base s)."""
    m_val, adrs = multi[k, s]
    return sum(
        np.exp(1j * 2 * np.pi * np.dot(q_red, svecs[adrs + ll]))
        for ll in range(m_val)
    ) / m_val

def dynamical_matrix(q, B):
    """
    D(q) 6×6 para Si diamante.

    Estructura de bloques:
        D = [ D_AA  D_AB ]   filas 0:3 → atomo A desplazado, svecs ref=0
            [ D_BA  D_BB ]   filas 3:6 → atomo B desplazado, svecs ref=1

    Separacion por subred:
        m par   (0,2,4,...) → subred A → columnas 0:3
        m impar (1,3,5,...) → subred B → columnas 3:6
    """
    q_red = q
    D = np.zeros((6, 6), dtype=complex)

    for m in range(n_atoms):
        ph_A = phase_avg(q_red, svecs, multi, m, s=0)
        ph_B = phase_avg(q_red, svecs, multi, m, s=1)

        if m % 2 == 0:   # atomo m es subred A → columnas 0:3
            D[0:3, 0:3] += IFC_A[m] * ph_A   # D_AA
            D[3:6, 0:3] += IFC_B[m] * ph_B   # D_BA
        else:             # atomo m es subred B → columnas 3:6
            D[0:3, 3:6] += IFC_A[m] * ph_A   # D_AB
            D[3:6, 3:6] += IFC_B[m] * ph_B   # D_BB

    D = (D + D.conj().T) / 2.0
    return D / M_Si

def phonon_frequencies(D):
    w2, _ = np.linalg.eigh(D)
    w2[w2 < 0] = 0.0
    return np.sqrt(w2) * factor

# ── Diagnóstico en Gamma ─────────────────────────────────────────
print("\n" + "=" * 60)
print("DIAGNOSTICO D(q) EN GAMMA q=[0,0,0]")
print("=" * 60)
D_gamma = dynamical_matrix(np.array([0.0, 0.0, 0.0]), B)
print("  Parte real:")
print("  " + str(np.round(D_gamma.real, 6)).replace("\n", "\n  "))
print(f"  Norma parte imaginaria (debe ser ~0): {np.linalg.norm(D_gamma.imag):.2e}")
freqs_gamma = phonon_frequencies(D_gamma)
print(f"  Frecuencias en Gamma (cm⁻¹): {np.round(freqs_gamma, 4)}")
print(f"  (Esperado: 3 modos acusticos ~0, 3 opticos ~520 cm⁻¹)")

# ── Diagnóstico en X ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("DIAGNOSTICO D(q) EN X q=[0.5,0.5,0]")
print("=" * 60)
D_X = dynamical_matrix(np.array([0.5, 0.5, 0.0]), B)
freqs_X = phonon_frequencies(D_X)
print(f"  Frecuencias en X (cm⁻¹): {np.round(freqs_X, 4)}")
print(f"  (Referencia DFPT Si en X: ~150, ~150, ~400, ~400, ~430, ~430 cm⁻¹)")

# ── Diagnóstico en L ─────────────────────────────────────────────
print("\n" + "=" * 60)
print("DIAGNOSTICO D(q) EN L q=[0.5,0.5,0.5]")
print("=" * 60)
D_L = dynamical_matrix(np.array([0.5, 0.5, 0.5]), B)
freqs_L = phonon_frequencies(D_L)
print(f"  Frecuencias en L (cm⁻¹): {np.round(freqs_L, 4)}")

# ═════════════════════════════════════════════════════════════════
#  CALCULAR BANDAS
# ═════════════════════════════════════════════════════════════════

print("\n" + "=" * 60)
print("CALCULANDO BANDAS...")
print("=" * 60)

omega = np.zeros((QN, n_branch))
for i in range(QN):
    omega[i] = phonon_frequencies(dynamical_matrix(q_path_red[i], B))

n_neg = np.sum(np.any(omega == 0.0, axis=1))
print(f"  Puntos q con algun modo forzado a 0: {n_neg}")
print(f"  Frecuencia maxima : {omega.max():.4f} cm⁻¹")
print(f"  Frecuencia minima no nula: {omega[omega > 1.0].min():.4f} cm⁻¹")

bandas = np.column_stack((s_path, omega))
# ═════════════════════════════════════════════════════════════════
#  GRAFICA
# ═════════════════════════════════════════════════════════════════

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
plt.ylabel(r'$\omega$ (cm⁻¹)')
plt.xticks(tick_positions, tick_labels)
plt.tight_layout()
plt.savefig(output_plot)
print(f"\n  Grafica guardada en {output_plot}")

np.savetxt(output_dat, bandas)
print(f"  Datos guardados en {output_dat}")
