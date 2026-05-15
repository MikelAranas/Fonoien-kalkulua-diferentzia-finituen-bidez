import numpy as np
import re

def read_forces_qe(out_file, nat):
    """
    Lee las fuerzas finales de un output de pw.x
    Devuelve array (nat, 3) en Ry/Bohr
    """
    forces = np.zeros((nat, 3))
    with open(out_file, "r") as f:
        lines = f.readlines()

    for i, line in enumerate(lines):
        if "Forces acting on atoms" in line:
            for a in range(nat):
                l = lines[i + 2 + a]
                nums = re.findall(r"[-+]?\d*\.\d+E?[+-]?\d*", l)
                forces[a, :] = [float(nums[-3]),
                                float(nums[-2]),
                                float(nums[-1])]
            return forces

    raise RuntimeError(f"No se encontraron fuerzas en {out_file}")

nat = 27
forces_plus  = {} 
forces_minus = {}

directions = {"x": 0, "y": 1, "z": 2}
i=0 #bakarrik lehenengo atomoa desplazatuz
for mu in directions:
    n=i
    key = (i, directions[mu],n)  # (i=0, mu=desplazamenduaren norabidea,n=0), kasu honetan atomoagaz bat, gelaxka bakoitzak atomo bakarra dekolako bestela n=j/atomo kop.
    f_plus  = f"atom{i+1:02d}_{mu}_plus.out" 
    f_minus = f"atom{i+1:02d}_{mu}_minus.out"
    forces_plus[key]  = np.array(read_forces_qe(f_plus, nat))#hiztegizen (i-1,mu,n) giltzegaz gordetan deu, i garren atomoa mu norabidean desplazatzerakoan atomoek jasandako indarren lista bat
    forces_minus[key] = np.array(read_forces_qe(f_minus, nat))#bardin, baina aurkako noranzkoan desplazatzean
#HEMETIK hiru elemntuko bi hiztegi

DELTA=0.01889726126 #hau bohr-etan erabili dan desplazamendue

Phi = np.zeros((3, 81))  # Ry / Bohr^2

for mu in range(3):#0,1,2 noranzkoa marketan deu. (j=0,beta=0) lehenengo atomoak x norabidien jasan dabezen indarrakaz
    row = 3*i + mu
    n=i
    Fp = forces_plus[(i, mu, n)] 
    Fm = forces_minus[(i, mu, n)]
    dF = (Fp - Fm) / (2 * DELTA)  # Ry / Bohr^2 # IFC-n filie
    for j in range(nat): # Oin dF filako elementuek banaka kolokeu j,nu indizien arabera
        for nu in range(3):
            col = 3*j + nu
            Phi[row, col] = -dF[j, nu]

for mu in range(3):
    for nu in range(3):
        col_sum = sum(Phi[mu, 3*j + nu] for j in range(nat))
        correction = col_sum / nat
        for j in range(nat):
            Phi[mu, 3*j + nu] -= correction

# APLICAR ACOUSTIC SUM RULE
#for mu in range(3):
#    row = 3*0 + mu
#    # La suma de todas las columnas en esta fila debe ser 0
#    row_sum = np.sum(Phi[row, :])
#    # Distribuir la corrección sobre todas las columnas
#    correction = row_sum / 24
#    Phi[row, :] -= correction

np.savetxt("IFC_murriztue.dat", Phi)
np.set_printoptions(suppress=True, precision=8)
print(Phi)
