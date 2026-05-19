import matplotlib as mpl
import matplotlib.pyplot as plt
import numpy as np
from numpy import linalg

# Cargar frecuencias (columna 0 de QE la ignoramos, usamos nuestro s_path)
freq = np.loadtxt('si.freq.gp', unpack=True, usecols=(0,1,2,3,4,5,6))

# Cargar path en coordenadas cristalinas
q_path = np.loadtxt('path.dat', unpack=False, usecols=(0,1,2))

# Vectores de la red reciproca del FCC (en unidades de 2pi/a)
b1 = np.array([-1.0, -1.0,  1.0])
b2 = np.array([ 1.0,  1.0,  1.0])
b3 = np.array([-1.0,  1.0, -1.0])
B = 0.707107*np.vstack([b1, b2, b3])

# Convertir path a cartesianas y calcular distancia acumulada
QN = q_path.shape[0]
q_path_cart = q_path @ B
s_path = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_path_cart[i] - q_path_cart[i-1])
# Guardar bandak_cm.dat: primera columna s_path, resto frecuencias de si.freq.gp
# freq[0] es el eje x original, freq[1:] son las bandas
# Sustituimos freq[0] por s_path
freq[1:]=freq[1:]/33.3564
data_out = np.column_stack([freq[i] for i in range(len(freq))])
np.savetxt('bandak.dat', data_out, fmt='%12.6f')
# Define high-symmetry points
G=[0.0,0.0,0.0]
X=[0.5,0.5,0.0]
K=[0.375,0.625,0.0]
G2=[0.0,1.0,0.0]
L=[0.0,1.5,0.0]
sympoints=[G,X,K,G2,L]

# Start figure
fig, ax = plt.subplots(figsize=(8*(1+np.sqrt(5)/2),10*1))
# Plot phonon bands
for i in range(1,len(freq)):
    ax.plot(freq[0], freq[i], color='blue')
# Set maximum energy of plot
wmax = np.max(freq)
ytop = wmax + 2.0

# Set minimum energy of plot
wmin = np.min(freq)
ybot = wmin - 0.01
# If negative modes, highlight region
if (wmin < 0.0):
    ax.plot([np.min(freq[0]),np.max(freq[0])],[0.0,0.0], color='black', alpha=0.5)
    ax.fill_between(freq[0], 0.0, ybot, color='gray', alpha=0.5)

# Plot high-symmetry points
for i in range(0,len(q_path)):
    for q in sympoints:
        if (np.linalg.norm(q_path[i]-q)<1E-05):
            ax.plot([freq[0][i],freq[0][i]],[ybot,ytop],'--',color='gray',linewidth=0.7)
           # Name the symmetry point
            if (q == G):
                ax.text(freq[0][i],ybot-25.0,r'$\Gamma$',fontsize=25)
            if (q == X):
                ax.text(freq[0][i],ybot-25.0,r'$\mathrm{X}$',fontsize=25)
            if (q == K):
                ax.text(freq[0][i],ybot-25.0,r'$\mathrm{K}$',fontsize=25)
            if (q == G2):
                ax.text(freq[0][i],ybot-25.0,r'$\Gamma$',fontsize=25)
            if (q == L):
                ax.text(freq[0][i],ybot-25.0,r'$\mathrm{L}$',fontsize=25)

# Write label on y axis
ax.set_ylabel(r'$\omega \/ (\mathrm{THz})$',fontsize=30,labelpad=10)

# Set plot limits
ax.set_xlim(np.min(freq[0]),np.max(freq[0]))
ax.set_ylim(ybot,ytop)


# Set tick parameters
ax.tick_params(axis='x', color='black', labelsize='0', pad=0, length=0, width=0)
ax.tick_params(axis='y', color='black', labelsize='25', pad=5, length=5, width=2)

plt.savefig('DFPT_4x4x4.pdf')
