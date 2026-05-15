import numpy as np
from numpy import linalg
from collections import defaultdict
import cmath
import matplotlib.pyplot as plt
import sys, os
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                 '..', '..', '..', '..', 'modules'))
from svecs_module import get_svecs_multi



# q_path-a kargatzen da aukertautako fitxategitik (printzipioz DFPT-n erabili den ibilbide mantentzeko)
# Ohartu: koordenatu kristalinotan ari garela lan egiten
q_path = np.loadtxt('path.dat', skiprows=0, usecols=(0,1,2))
QN = q_path.shape[0]


# Koordenatu kartestarretara aldatzeko:
# pw.x-eko output fitxategitik elkarrekiko sareko oinarria atera
# Ohartu: 2pi/alat unitatetan dauedela, oinarri kartestarrean kx,ky,kz
b1 = np.array([-0.707107, -0.707107,  0.707107])
b2 = np.array([ 0.707107,  0.707107,  0.707107])
b3 = np.array([-0.707107,  0.707107, -0.707107])

# B matrizea, b oinarri bektoreak lerroetan
B = np.vstack([b1, b2, b3])  # shape (3,3)
# q_cart = q1*b1 + q2*b2 + q3*b3
q_path2 = q_path @ B   # shape (QN,3)


# Orain, s_path, ibilbidean akumulatzen den distantzia puntu bakoitzerako kalkulatuko dugu
# Ohartu: script honetan koordenatu kartestarretan lan egingo dugula
s_path = np.zeros(QN)
for i in range(1, QN):
    s_path[i] = s_path[i-1] + linalg.norm(q_path2[i] - q_path2[i-1])

# Behin x ardatza prest dagoela
# y ardatza kalkulatuko da, q bakoitzari dagozkion frekuentziak hain zuzen

#Horretarako lehenengo IFC murriztua kargatu
IFC = np.loadtxt("IFC_murriztue.dat")

# Atomo bakoitza aurkitzen den gelaxkari dagokion bektoreen numpy array bat (faserako)
# Ohartu: koordenatu kristalinotan idatzi ditugula orain

# Koordenatu kartestarretara aldatzeko
# Horretarako .out fitxategitik sare zuzeneko oinarri bektoreak erabiliko ditugu, hauek alatetan egonik, eta oinarri kartestarrean x,y,z
a1 = ([-0.707107,0.000000,0.707107])
a2 = ([0.000000,0.707107,0.707107])
a3 = ([-0.707107,0.707107,0.000000])

A =np.vstack([a1,a2,a3])

N = 3
# sc_positions in primitive fractional coords, in the same order as QE
# ATOMIC_POSITIONS (ix outermost, iz innermost).
sc_positions = np.array(
    [[ix, iy, iz] for ix in range(N) for iy in range(N) for iz in range(N)],
    dtype=float,
)
svecs, multi = get_svecs_multi(A, sc_positions, N)

IFC_m = defaultdict(lambda: np.zeros((3,3),dtype=float))# 0z beteriko 3x3 matrizeak dituen hiztegia.

#{0:(3,3),1:(3,3),...} bakarrik behar ditugunak  matrize dinamikoaren batukarirako
n=0
for m in range(27):
    for mu in range(3):
        for nu in range(3):
            row = 3*n + mu #n=0 dalez bakarrik n=0 bloke filan geratzen da batukarize, IFC Matrizien lehenengo hiru errenkadak eta 24 zutabiek, 8 matrizetan banatute.
            col = 3*m  + nu
            IFC_m[m][mu, nu] = IFC[row, col]#*778.35292 #R_list-en batutako posizioei dagozkien matrizeak eraiki. Horretarako IFC-ko lehen hiru errenkadak eta zutabe danak erabiliz

M_Al = 26.9815
# Unidades de la matriz dinámica: IFC en Ry/bohr², masa en uma
# -> autovalores w² en Ry/(bohr²·uma)
#
# Paso 1: Ry/bohr² -> eV/Å²
#   1 Ry = 13.605693 eV
#   1 bohr = 0.529177 Å  ->  1 bohr² = 0.529177² Å²
#   factor parcial: sqrt(13.605693 / 0.529177²)
#
# Paso 2: eV/Å²/uma -> THz
#   1 eV = 1.60218e-19 J
#   1 uma = 1.66054e-27 kg
#   1 Å = 1e-10 m
#   sqrt(eV/Å²/uma) = sqrt(J/m²/kg) = sqrt(s⁻²) = rad/s
#   dividiendo por 2π*1e12 se obtienen THz
#   -> el factor numérico resultante es 15.633302
#
# Paso 3: THz -> cm⁻¹
#   ṽ = ν/c,  c = 2.99792458e10 cm/s
#   1 THz = 1e12 Hz / 2.99792458e10 cm/s = 33.3564 cm⁻¹
#
# Factor completo: Ry/bohr²/uma -> eV/Å²/uma -> THz -> cm⁻¹
factor = np.sqrt(13.605693 / 0.529177**2) * 15.633302 * 33.3564

def dynamical_matrix(q, IFC_m):
    """
    q : vector q en coordenadas reducidas
    IFC_m : dict {m: 3x3 IFC}
    """
    q_red = q
    
    D = np.zeros((3,3), dtype=complex)
    
    # Loop sobre átomos de supercelda
    for k in range(27):
        m_val, adrs = multi[k, 0]  # Multiplicidad
        
        # Promediar fase sobre vectores equivalentes
        phase_sum = 0.0
        for ll in range(m_val):
            vec = svecs[adrs + ll]
            phase_sum += np.exp(1j * 2 * np.pi * np.dot(q_red, vec))
        
        phase_avg = phase_sum / m_val
        D += IFC_m[k] * phase_avg
    
    # Hermitianizar
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


###############
#   GRAFIKA   #
###############

# Funtzio hau erabiliz, Labelak spath puntu egokizetan ahalkodire jarri
def s_at_q_on_polyline(q_path, s_path, q_target):
    """
    Devuelve la posición 's' a lo largo del polilínea (q_path, s_path)
    donde el punto q_target proyecta más cerca. Preciso incluso si q_target
    no coincide con ningún nodo de q_path.
    """
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

# Γ–X–K–Γ*–L
# Ohartu, koordenatu kristalinotan gaudela!, 2pi/alatetan
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

# Puntu berezietan lerro bertikalak
for s_tick in tick_positions:
    plt.axvline(s_tick, color='0.7', ls=':', lw=0.8, zorder=0)

plt.xlim(s_path[0], s_path[-1])
plt.ylabel(r'$\omega$ (cm⁻¹)')
plt.xticks(tick_positions, tick_labels)
plt.tight_layout()
plt.savefig("Al_DF_bandak.pdf", dpi=200)

# bandak.dat-en gordeko dire bai x ardatzeko balioak bai y, DFPT-ga komparaketa eitzeko
np.savetxt("Al_DF_bandak.dat",bandak)




