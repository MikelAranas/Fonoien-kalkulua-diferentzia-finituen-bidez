import numpy as np
from numpy import linalg
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', '..', '..', 'modules'))
from svecs_module import get_svecs_multi

# ═══════════════════════════════════════════════════════════════
#  PARAMETROS
# ═══════════════════════════════════════════════════════════════

M_C     = 12.011
N_BASIS = 2
N       = 2

IFC_FILE  = 'IFC_murriztue.dat'
PATH_FILE = 'path.dat'
OUTPUT_PDF = 'Grafenoa_DF_bandak.pdf'
OUTPUT_DAT = 'Grafenoa_DF_bandak.dat'

# Vectores primitivos (del c.scf.out, en alat)
a1 = np.array([ 1.000000,  0.000000,  0.000000])
a2 = np.array([-0.500000,  np.sqrt(3)/2,  0.000000])
a3 = np.array([ 0.000000,  0.000000,  3.610072])
A_PRIM = np.vstack([a1, a2, a3])

# Vectores reciprocos (del c.scf.out, en 2pi/alat)
b1 = np.array([1.000000,  0.577350, 0.000000])
b2 = np.array([0.000000,  1.154701, 0.000000])
b3 = np.array([0.000000,  0.000000, 0.277003])
B = np.vstack([b1, b2, b3])

factor = np.sqrt(13.605693 / 0.529177**2) * 15.633302 * 33.3564

symmetry_points = [
    (r'$\Gamma$', np.array([0.000000, 0.000000, 0.0])),
    (r'$M$',      np.array([0.500000, 0.000000, 0.0])),
    (r'$K$',      np.array([0.333333, 0.333333, 0.0])),
    (r'$\Gamma$', np.array([0.000000, 0.000000, 0.0])),
]

# ═══════════════════════════════════════════════════════════════
#  POSICIONES EN COORDENADAS PRIMITIVAS (2x2x1, 8 atomos)
# ═══════════════════════════════════════════════════════════════

sc_positions = np.array([
    [0,      0,      0.0],   # 0  subred A
    [1/3,    2/3,    0.0],   # 1  subred B
    [1,      0,      0.0],   # 2  A
    [1+1/3,  2/3,    0.0],   # 3  B  — exact 4/3 for correct mult=3 in svecs
    [0,      1,      0.0],   # 4  A
    [1/3,    1+2/3,  0.0],   # 5  B  — exact 5/3
    [1,      1,      0.0],   # 6  A
    [1+1/3,  1+2/3,  0.0],   # 7  B  — exact 4/3, 5/3
], dtype=float)

nat = len(sc_positions)
print(f"N={N}, nat={nat}")

# ═══════════════════════════════════════════════════════════════
#  SVECS Y MULTI
# ═══════════════════════════════════════════════════════════════

basis_positions = np.array([
    [0.0,    0.0,    0.0],   # C_A
    [1/3,    2/3,    0.0],   # C_B
])
svecs, multi = get_svecs_multi(A_PRIM, sc_positions, [N, N, 1], basis_positions=basis_positions)

# ═══════════════════════════════════════════════════════════════
#  CARGA DE IFC
# ═══════════════════════════════════════════════════════════════

IFC = np.loadtxt(IFC_FILE)
print(f"IFC shape: {IFC.shape}  (esperado: ({3*N_BASIS}, {3*nat}))")

IFC_blocks = {}
for alpha in range(N_BASIS):
    IFC_blocks[alpha] = {}
    for m in range(nat):
        IFC_blocks[alpha][m] = IFC[3*alpha:3*alpha+3, 3*m:3*m+3]

# ═══════════════════════════════════════════════════════════════
#  Q-PATH
# ═══════════════════════════════════════════════════════════════

q_path = np.loadtxt(PATH_FILE, usecols=(0,1,2))
QN = q_path.shape[0]
q_cart = q_path @ B
s_path = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_cart[i] - q_cart[i-1])

# ═══════════════════════════════════════════════════════════════
#  MATRIZ DINAMICA
# ═══════════════════════════════════════════════════════════════

def dynamical_matrix(q_red):
    D = np.zeros((3*N_BASIS, 3*N_BASIS), dtype=complex)
    for alpha in range(N_BASIS):
        for m in range(nat):
            beta = m % N_BASIS
            m_val, adrs = multi[m, alpha]
            phase = sum(
                np.exp(1j * 2 * np.pi * np.dot(q_red, svecs[adrs + ll]))
                for ll in range(m_val)
            ) / m_val
            D[3*alpha:3*alpha+3, 3*beta:3*beta+3] += IFC_blocks[alpha][m] * phase
    D = (D + D.conj().T) / 2.0
    return D / M_C

def phonon_frequencies(D):
    w2, _ = np.linalg.eigh(D)
    w2[w2 < 0] = 0.0
    return np.sqrt(w2) * factor

# Diagnostico en Gamma
print("\nFrecuencias en Gamma:")
freqs_G = phonon_frequencies(dynamical_matrix(np.array([0.0, 0.0, 0.0])))
print(f"  {np.round(freqs_G, 2)} cm-1")
print("  (Esperado: 3 acusticos ~0, ZO ~900, LO/TO ~1580 cm-1)")

print("\nFrecuencias en M:")
freqs_M = phonon_frequencies(dynamical_matrix(np.array([0.5, 0.0, 0.0])))
print(f"  {np.round(freqs_M, 2)} cm-1")

# ═══════════════════════════════════════════════════════════════
#  CALCULO DE BANDAS
# ═══════════════════════════════════════════════════════════════

n_branch = 3 * N_BASIS
omega = np.zeros((QN, n_branch))
for i in range(QN):
    omega[i] = phonon_frequencies(dynamical_matrix(q_path[i]))

print(f"Frecuencia maxima: {omega.max():.2f} cm-1")

# ═══════════════════════════════════════════════════════════════
#  GRAFICA
# ═══════════════════════════════════════════════════════════════

def s_at_q(q_path, s_path, q_target):
    best_s, best_d2 = None, np.inf
    for i in range(len(q_path)-1):
        a, b = q_path[i], q_path[i+1]
        v = b - a
        vv = np.dot(v, v)
        if vv == 0: continue
        t = np.clip(np.dot(q_target - a, v) / vv, 0, 1)
        d2 = np.dot(q_target - (a + t*v), q_target - (a + t*v))
        if d2 < best_d2:
            best_d2 = d2
            best_s = s_path[i] + t * linalg.norm(v)
    return best_s

tick_labels, tick_positions = [], []
for label, qpt in symmetry_points:
    tick_labels.append(label)
    tick_positions.append(s_at_q(q_cart, s_path, qpt @ B))

fig, ax = plt.subplots(figsize=(7.5, 4.5))
for b in range(n_branch):
    ax.plot(s_path, omega[:, b], color='C0', lw=1.5)
for s_tick in tick_positions:
    ax.axvline(s_tick, color='0.7', ls=':', lw=0.8)
ax.set_xlim(s_path[0], s_path[-1])
ax.set_ylim(bottom=0)
ax.set_ylabel(r'$\omega$ (cm$^{-1}$)')
ax.set_xticks(tick_positions)
ax.set_xticklabels(tick_labels)
plt.tight_layout()
plt.savefig(OUTPUT_PDF)
print(f"Grafica guardada: {OUTPUT_PDF}")

np.savetxt(OUTPUT_DAT, np.column_stack((s_path, omega)))
print(f"Datos guardados: {OUTPUT_DAT}")

