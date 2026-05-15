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

nat = 432
forces_plus  = {} 
forces_minus = {}

directions = {"x": 0, "y": 1, "z": 2}

for i in range(2):
    for mu in directions:
        n=0
        key = (i, directions[mu],n)  # (i=0, mu=desplazamenduaren norabidea,n=0, kasu honetan atomoagaz bat, gelaxka bakoitzak atomo bakarra dekolako bestela n=j/atomo kop.)
        f_plus  = f"atom{i+1:02d}_{mu}_plus.out" 
        f_minus = f"atom{i+1:02d}_{mu}_minus.out"
        forces_plus[key]  = np.array(read_forces_qe(f_plus, nat))#hiztegizen (i-1,mu) giltzegaz gordetan deu, i garren atomoa mu norabidean desplazatzerakoan atomoek jasandako indarren lista bat
        forces_minus[key] = np.array(read_forces_qe(f_minus, nat))#bardin, baina aurkako noranzkoan desplazatzean
#HEMETIK hiru elemntuko bi hiztegi
DELTA_ang  = 0.01                        # Å, el que usaste en crearinputs.py
DELTA= DELTA_ang / 0.529177        # Bohr, para dividir las fuerzas de QE

Phi = np.zeros((6, 1296))  # Ry / Bohr^2
for i in range(2):
    for mu in range(3):#0,1,2 noranzkoa marketan deu. (j=0,beta=0) lehenengo atomoak x norabidien jasan dabezen indarrakaz
        row = 3*i + mu
        Fp = forces_plus[(i, mu, n)] 
        Fm = forces_minus[(i, mu, n)]
        dF = (Fp - Fm) / (2 * DELTA)  # Ry / Bohr^2 # IFC-n filie
        for j in range(nat): # Oin dF filako elementuek banaka kolokeu j,nu indizien arabera
            for nu in range(3):
                col = 3*j + nu
                Phi[row, col] = -dF[j, nu]
for i in range(2):
    for mu in range(3):
        for nu in range(3):
            row = 3*i + mu
            col_sum = sum(Phi[row, 3*j + nu] for j in range(nat))
            correction = col_sum / nat
            for j in range(nat):
                Phi[row, 3*j + nu] -= correction

print("Sumas de fila despues de ASR (deben ser ~0):")
for row in range(6):
    print(f"  fila {row}: {np.sum(Phi[row,:]):.2e}")

np.savetxt("IFC_murriztue.dat", Phi)
