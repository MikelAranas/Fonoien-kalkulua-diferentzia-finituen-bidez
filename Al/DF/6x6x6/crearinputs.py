#!/usr/bin/env python3
import os
import re
import numpy as np
from os.path import expanduser, abspath

# ==========================
# Parámetros
# ==========================
# Ajusta si tu base está en otra ruta
INPUT_BASE = './al_sup.scf.in'
OUTPUT_DIR = './in'
DELTA = 0.01  # Å
DIRECTIONS = ['x', 'y', 'z']  # cartesianas

# Nombre del prefix (debe coincidir con el SCF base)
PREFIX = 'al'

# Resolver ruta del input base
INPUT_BASE = abspath(expanduser(INPUT_BASE))
os.makedirs(OUTPUT_DIR, exist_ok=True)

print(f"[INFO] Leyendo input base: {INPUT_BASE}")
if not os.path.isfile(INPUT_BASE):
    raise FileNotFoundError(f"No se encuentra el archivo base: {INPUT_BASE}")

# ==========================
# Leer archivo base
# ==========================
with open(INPUT_BASE, 'r') as f:
    text = f.read()

# Extraer CELL_PARAMETERS (angstrom)
cell_match = re.search(r'CELL_PARAMETERS\s+angstrom\s*(.*?)\n\s*\n', text, re.S)
if not cell_match:
    cell_match = re.search(r'CELL_PARAMETERS\s+angstrom\s*(.*?)(?:\n[A-Z_]+\b)', text, re.S)
if not cell_match:
    raise ValueError("No se encontró bloque CELL_PARAMETERS angstrom en el input base.")

cell_block = [l for l in cell_match.group(1).strip().splitlines() if l.strip()]
if len(cell_block) < 3:
    raise ValueError("Bloque CELL_PARAMETERS incompleto (menos de 3 líneas).")

A = np.array([[float(x) for x in re.split(r'\s+', line.strip())] for line in cell_block[:3]])
print("[INFO] Matriz de celda (A):")
print(A)

# Detectar bloque ATOMIC_POSITIONS y su unidad
apos_match = re.search(r'ATOMIC_POSITIONS\s+(\w+)\s*(.*?)\n\s*\n', text, re.S)
if not apos_match:
    apos_match = re.search(r'ATOMIC_POSITIONS\s+(\w+)\s*(.*?)(?:\nK_POINTS\b|\nCELL_PARAMETERS\b|\n[A-Z_]+\b)', text, re.S)
if not apos_match:
    raise ValueError("No se encontró bloque ATOMIC_POSITIONS en el input base.")

apos_unit = apos_match.group(1).lower()  # 'crystal' o 'angstrom'
apos_lines = [l for l in apos_match.group(2).strip().splitlines() if l.strip()]

atoms, coords = [], []
for line in apos_lines:
    parts = re.split(r'\s+', line.strip())
    if len(parts) < 4:  # ignora líneas vacías/comentarios
        continue
    atom = parts[0]
    xyz = np.array(list(map(float, parts[1:4])))
    atoms.append(atom)
    coords.append(xyz)
coords = np.vstack(coords)
N = len(atoms)
print(f"[INFO] ATOMIC_POSITIONS unidad: {apos_unit} | átomos detectados: {N}")

# Convertir a cartesianas si están en 'crystal'
if apos_unit == 'crystal':
    # r_cart = frac.dot(A)  (A con vectores por filas)
    R_cart = coords.dot(A)
else:
    R_cart = coords.copy()

# ==========================
# Generar inputs desplazados (forzando ATOMIC_POSITIONS angstrom)
# ==========================
axis_map = {'x': 0, 'y': 1, 'z': 2}
count = 0
ATOMS_TO_DISPLACE = [0]  # En el aluminio solo desplazaremos el primer atomo

for i_atom in ATOMS_TO_DISPLACE:
    for dir_name in DIRECTIONS:
        a = axis_map[dir_name]
        for sign, sname in [(1, 'plus'), (-1, 'minus')]:
            disp = np.zeros(3)
            disp[a] = sign * DELTA
            R_new = R_cart.copy()
            R_new[i_atom, :] = R_cart[i_atom, :] + disp

            # Construir bloque nuevo en angstrom
            apos_block_new = ["ATOMIC_POSITIONS angstrom"]
            for k in range(N):
                apos_block_new.append(f"{atoms[k]} {R_new[k,0]:.6f} {R_new[k,1]:.6f} {R_new[k,2]:.6f}")
            apos_block_new_str = "\n".join(apos_block_new) + "\n\n"

            # Reemplazar ATOMIC_POSITIONS por el nuevo bloque
            new_text = re.sub(r'ATOMIC_POSITIONS\s+\w+\s*.*?\n\s*\n',
                              apos_block_new_str, text, flags=re.S)
            new_text = re.sub(r"prefix\s*=\s*'.*?'",     f"prefix = '{PREFIX}'",        new_text)

            fname = f"atom{i_atom+1:02d}_{dir_name}_{sname}.in"
            outpath = os.path.join(OUTPUT_DIR, fname)
            with open(outpath, 'w') as fout:
                fout.write(new_text)
            count += 1

print(f"[OK] Generados {count} archivos en '{abspath(OUTPUT_DIR)}'")
print("[Siguiente] Usa el launcher para ejecutar en paralelo con scratch por job y evitar colisiones.")

