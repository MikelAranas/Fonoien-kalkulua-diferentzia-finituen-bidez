import os
import shutil

HERE   = os.path.dirname(os.path.abspath(__file__))
MAT    = os.path.dirname(HERE)
OUT_DIR  = HERE
DF_DIR   = os.path.join(MAT, 'DF')
DFPT_DIR = os.path.join(MAT, 'DFPT')

sizes = sorted(d for d in os.listdir(DF_DIR)
               if os.path.isdir(os.path.join(DF_DIR, d)))
print(f"Supercell sizes found: {sizes}")

mat_name = os.path.basename(MAT)
df_band_name = f'{mat_name}_DF_bandak.dat'

dfpt_band_candidates = ['bandak.dat', 'graphene_bandak_cm.dat']

for size in sizes:
    src = os.path.join(DF_DIR, size, 'out', df_band_name)
    dst = os.path.join(OUT_DIR, f'DF_{size}.dat')
    if os.path.isfile(src):
        shutil.copy2(src, dst)
        print(f"  DF   {size}: {src} -> {dst}")
    else:
        print(f"  DF   {size}: missing {src}")
    for name in dfpt_band_candidates:
        src = os.path.join(DFPT_DIR, size, name)
        if os.path.isfile(src):
            break
    else:
        print(f"  DFPT {size}: missing (tried {dfpt_band_candidates})")
        continue
    dst = os.path.join(OUT_DIR, f'DFPT_{size}.dat')
    shutil.copy2(src, dst)
    print(f"  DFPT {size}: {src} -> {dst}")
