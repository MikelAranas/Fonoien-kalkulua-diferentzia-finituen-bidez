import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg

freq = np.loadtxt('graphene.freq.gp', unpack=True)

ca = 8.0
b1 = np.array([1.000000,  0.577350, 0.000000])
b2 = np.array([0.000000,  1.154701, 0.000000])
b3 = np.array([0.000000,  0.000000, 0.232558])
B  = np.vstack([b1, b2, b3])

q_path = np.loadtxt('path.dat', unpack=False, usecols=(0, 1, 2))
QN = q_path.shape[0]

q_cart = q_path @ B
s_path = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_cart[i] - q_cart[i-1])
print(s_path)
data_out = np.column_stack([s_path] + [freq[i] for i in range(1, len(freq))])
np.savetxt('graphene_bandak_cm.dat', data_out, fmt='%12.6f')
print("Guardado graphene_bandak_cm.dat")

G = np.array([ 0.000000,  0.000000, 0.0])
M = np.array([ 0.000000, -0.577588, 0.0])
K = np.array([-0.333471, -0.577588, 0.0])
sympoints = [G, M, K]
symlabels = {
    tuple(G): r'$\Gamma$',
    tuple(M): r'$\mathrm{M}$',
    tuple(K): r'$\mathrm{K}$',
}

fig, ax = plt.subplots(figsize=(8 * (1 + np.sqrt(5) / 2), 10))

for i in range(1, len(freq)):
    ax.plot(s_path, freq[i], color='blue', lw=1.2)

wmax = np.max(freq[1:])
wmin = np.min(freq[1:])
ytop = wmax + 2.0
ybot = wmin - 0.01

if wmin < 0.0:
    ax.plot([s_path[0], s_path[-1]], [0.0, 0.0], color='black', alpha=0.5)
    ax.fill_between(s_path, 0.0, ybot, color='gray', alpha=0.5)

for i in range(QN):
    for q in sympoints:
        if linalg.norm(q_path[i] - q) < 1E-05:
            ax.plot([s_path[i], s_path[i]], [ybot, ytop], '--', color='gray', linewidth=0.7)
            ax.text(s_path[i], ybot - 25.0, symlabels[tuple(q)], fontsize=25, ha='center')

ax.set_ylabel(r'$\omega\ (\mathrm{cm^{-1}})$', fontsize=30, labelpad=10)
ax.set_xlim(s_path[0], s_path[-1])
ax.set_ylim(ybot, ytop)
ax.tick_params(axis='x', labelsize=0, length=0, width=0)
ax.tick_params(axis='y', labelsize=25, pad=5, length=5, width=2)

plt.tight_layout()
plt.savefig('graphene_bandak.pdf')
print("Guardado graphene_bandak.pdf")
