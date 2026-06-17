import numpy as np
from numpy import linalg
from collections import defaultdict
import cmath
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', '..', '..', 'modules'))
from svecs_module import get_svecs_multi

q_path = np.loadtxt('path.dat', skiprows=0, usecols=(0,1,2))
QN = q_path.shape[0]

b1 = np.array([-0.707107, -0.707107,  0.707107])
b2 = np.array([ 0.707107,  0.707107,  0.707107])
b3 = np.array([-0.707107,  0.707107, -0.707107])

B = np.vstack([b1, b2, b3])
q_path2 = q_path @ B

s_path = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_path2[i] - q_path2[i-1])

IFC = np.loadtxt("IFC_murriztue.dat")

a1 = ([-0.707107,0.000000,0.707107])
a2 = ([0.000000,0.707107,0.707107])
a3 = ([-0.707107,0.707107,0.000000])

A =np.vstack([a1,a2,a3])

N = 4
sc_positions = np.array(
    [[ix, iy, iz] for ix in range(N) for iy in range(N) for iz in range(N)],
    dtype=float,
)
svecs, multi = get_svecs_multi(A, sc_positions, N)

IFC_m = defaultdict(lambda: np.zeros((3,3),dtype=float))

n=0
for m in range(64):
    for mu in range(3):
        for nu in range(3):
            row = 3*n + mu
            col = 3*m  + nu
            IFC_m[m][mu, nu] = IFC[row, col]

M_Al = 26.9815
factor = np.sqrt(13.605693 / 0.529177**2) * 15.633302 * 33.3564

def dynamical_matrix(q, IFC_m):
    q_red = q
    
    D = np.zeros((3,3), dtype=complex)
    
    for k in range(64):
        m_val, adrs = multi[k, 0]
        
        phase_sum = 0.0
        for ll in range(m_val):
            vec = svecs[adrs + ll]
            phase_sum += np.exp(1j * 2 * np.pi * np.dot(q_red, vec))
        
        phase_avg = phase_sum / m_val
        D += IFC_m[k] * phase_avg
    
    D = (D + D.conj().T) / 2.0
    return D / M_Al

def phonon_frequencies(D):
    w2, _ = np.linalg.eigh(D)
    w2[w2 < 0] = 0.0
    return np.sqrt(w2) * factor

omega = np.zeros((QN, 3))
for i in range(QN):
    D = dynamical_matrix(q_path[i], IFC_m)
    omega[i] = phonon_frequencies(D)

bandak = np.column_stack((s_path, omega))

def s_at_q_on_polyline(q_path, s_path, q_target):
    best_s = None
    best_d2 = np.inf
    for i in range(len(q_path) - 1):
        a = q_path[i]
        b = q_path[i+1]
        v = b - a
        vv = np.dot(v, v)
        if vv == 0.0:
            continue
        t = np.dot(q_target - a, v) / vv
        t = np.clip(t, 0.0, 1.0)
        proj = a + t * v
        d2 = np.dot(q_target - proj, q_target - proj)
        if d2 < best_d2:
            best_d2 = d2
            seg_len = linalg.norm(v)
            best_s = s_path[i] + t * seg_len
    return best_s

tick_labels = []
tick_positions = []

p1 = ([0.0,   0.0, 0.0])
p2 = ([0.5,   0.5, 0.0])
p3 = ([0.375,0.625,0.0])
p4 = ([0.0,   1.0, 0.0])
p5 = ([0.0,   1.5, 0.0])

PATH = np.vstack([p1,p2,p3,p4,p5])

ordered = [(r'$\Gamma$',PATH[0]),(r'$X$',PATH[1]),(r'$K$',PATH[2]),(r'$\Gamma*$',PATH[3]),(r'$L$',PATH[4])]

for lab, qpt in ordered:
    s_star = s_at_q_on_polyline(q_path2, s_path, qpt @ B)
    tick_labels.append(lab)
    tick_positions.append(s_star)

plt.figure(figsize=(7.5, 4.0))
for b in range(3):
    plt.plot(s_path, omega[:, b], color='C0', lw=1.5)

for s_tick in tick_positions:
    plt.axvline(s_tick, color='0.7', ls=':', lw=0.8, zorder=0)

plt.xlim(s_path[0], s_path[-1])
plt.ylabel(r'$\omega$ (cm⁻¹)')
plt.xticks(tick_positions, tick_labels)
plt.tight_layout()
plt.savefig("Al_DF_bandak.pdf", dpi=200)

np.savetxt("Al_DF_bandak.dat",bandak)
