"""
Recoge los bandak.dat de cada supercell de DF/ y DFPT/ y los centraliza en
PLOTS/ con nombres claros.

Las rutas son relativas al directorio del script: detecta automáticamente
todas las supercells presentes (no hace falta listarlas a mano).
"""
import os
import shutil


HERE   = os.path.dirname(os.path.abspath(__file__))
MAT    = os.path.dirname(HERE)                # parent: Si/, Al/, Grafenoa/
PLOTS  = HERE
DF_DIR   = os.path.join(MAT, 'DF')
DFPT_DIR = os.path.join(MAT, 'DFPT')

# Auto-discover supercell sizes from DF/
sizes = sorted(d for d in os.listdir(DF_DIR)
               if os.path.isdir(os.path.join(DF_DIR, d)))
print(f"Supercell sizes found: {sizes}")

# DF bandak: in {MAT}/DF/{size}/out/{MAT_NAME}_DF_bandak.dat
mat_name = os.path.basename(MAT)               # Si | Al | Grafenoa
df_band_name = f'{mat_name}_DF_bandak.dat'

# DFPT bandak filename varies by material:
#   Si, Al  : bandak.dat                (in cm^-1)
#   Grafenoa: graphene_bandak_cm.dat    (in cm^-1, converted from QE matdyn .freq.gp)
dfpt_band_candidates = ['bandak.dat', 'graphene_bandak_cm.dat']

for size in sizes:
    # DF
    src = os.path.join(DF_DIR, size, 'out', df_band_name)
    dst = os.path.join(PLOTS, f'DF_{size}.dat')
    if os.path.isfile(src):
        shutil.copy2(src, dst)
        print(f"  DF   {size}: {src} -> {dst}")
    else:
        print(f"  DF   {size}: missing {src}")
    # DFPT
    for name in dfpt_band_candidates:
        src = os.path.join(DFPT_DIR, size, name)
        if os.path.isfile(src):
            break
    else:
        print(f"  DFPT {size}: missing (tried {dfpt_band_candidates})")
        continue
    dst = os.path.join(PLOTS, f'DFPT_{size}.dat')
    shutil.copy2(src, dst)
    print(f"  DFPT {size}: {src} -> {dst}")
