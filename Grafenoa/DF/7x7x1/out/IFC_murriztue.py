import numpy as np
import re

def read_forces_qe(out_file, nat):
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

nat = 98
n = 0
forces_plus  = {}
forces_minus = {}

directions = {"x": 0, "y": 1, "z": 2}

for i in range(2):
    for mu in directions:
        key = (i, directions[mu], n)
        f_plus  = f"atom{i+1:02d}_{mu}_plus.out"
        f_minus = f"atom{i+1:02d}_{mu}_minus.out"
        forces_plus[key]  = np.array(read_forces_qe(f_plus,  nat))
        forces_minus[key] = np.array(read_forces_qe(f_minus, nat))

DELTA_ang = 0.01
DELTA = DELTA_ang / 0.529177

Phi = np.zeros((6, 3 * nat))

for i in range(2):
    for mu in range(3):
        row = 3 * i + mu
        Fp = forces_plus[(i, mu, n)]
        Fm = forces_minus[(i, mu, n)]
        dF = (Fp - Fm) / (2 * DELTA)
        for j in range(nat):
            for nu in range(3):
                col = 3 * j + nu
                Phi[row, col] = -dF[j, nu]

for i in range(2):
    for mu in range(3):
        row = 3 * i + mu
        for nu in range(3):
            col_sum = sum(Phi[row, 3 * j + nu] for j in range(nat))
            correction = col_sum / nat
            for j in range(nat):
                Phi[row, 3 * j + nu] -= correction

print("Sumas de fila despues de ASR (deben ser ~0):")
for row in range(6):
    print(f"  fila {row}: {np.sum(Phi[row, :]):.2e}")

np.savetxt("IFC_murriztue.dat", Phi)
