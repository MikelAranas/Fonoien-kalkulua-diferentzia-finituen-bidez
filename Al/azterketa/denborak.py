#!/usr/bin/env python3
import os
import re
import csv

HERE = os.path.dirname(os.path.abspath(__file__))
HOME = os.path.dirname(HERE)
OUTPUT = os.path.join(HERE, 'denborak.csv')

calculos = {
    'DFPT': {
        '2x2x2': os.path.join(HOME, 'DFPT','2x2x2'),
        '3x3x3': os.path.join(HOME, 'DFPT','3x3x3'),
        '4x4x4': os.path.join(HOME, 'DFPT','4x4x4'),
        '6x6x6': os.path.join(HOME, 'DFPT','6x6x6'),
    },
    'DF': {
        '2x2x2': os.path.join(HOME, 'DF', '2x2x2'),
        '3x3x3': os.path.join(HOME, 'DF', '3x3x3'),
        '4x4x4': os.path.join(HOME, 'DF', '4x4x4'),
        '6x6x6': os.path.join(HOME, 'DF', '6x6x6')
    },
}

def parse_wall_time(time_str):
    total = 0.0
    d = re.search(r'(\d+\.?\d*)d', time_str)
    h = re.search(r'(\d+\.?\d*)h', time_str)
    m = re.search(r'(\d+\.?\d*)m', time_str)
    s = re.search(r'(\d+\.?\d*)s', time_str)
    if d: total += float(d.group(1)) * 86400
    if h: total += float(h.group(1)) * 3600
    if m: total += float(m.group(1)) * 60
    if s: total += float(s.group(1))
    return total

def find_wall_time(filepath):
    last = None
    try:
        with open(filepath, 'r', errors='ignore') as f:
            for line in f:
                m = re.search(r'(?:PWSCF|PHONON)\s*:\s*[\d\w\s.]+CPU\s+([\dhmsad. ]+)\s*WALL', line)
                if m:
                    last = parse_wall_time(m.group(1).strip())
    except:
        pass
    return last

def get_out_files(dirpath, method):
    out_files = []
    if method == 'DFPT':
        for fname in os.listdir(dirpath):
            if fname.endswith('.out') and ('ph' in fname or 'scf' in fname):
                out_files.append(os.path.join(dirpath, fname))
    elif method == 'DF':
        for fname in os.listdir(dirpath):
            if fname.endswith('.out') and 'scf' in fname.lower():
                out_files.append(os.path.join(dirpath, fname))
        out_dir = os.path.join(dirpath, 'out')
        if os.path.isdir(out_dir):
            for fname in sorted(os.listdir(out_dir)):
                if fname.endswith('.out'):
                    out_files.append(os.path.join(out_dir, fname))
    return out_files

rows = []
print(f"\n{'Method':<6} {'Size':<8} {'File':<45} {'WALL (s)':>10}")
print("-" * 75)

for method, sizes in calculos.items():
    for size, dirpath in sizes.items():
        if not os.path.isdir(dirpath):
            print(f"{method:<6} {size:<8} DIRECTORY NOT FOUND: {dirpath}")
            continue
        out_files = get_out_files(dirpath, method)
        total = 0.0
        for fpath in out_files:
            wall = find_wall_time(fpath)
            if wall is not None:
                fname_short = os.path.relpath(fpath, dirpath)
                print(f"{method:<6} {size:<8} {fname_short:<45} {wall:>10.1f}s")
                total += wall
        print(f"{method:<6} {size:<8} {'>>> TOTAL':<45} {total:>10.1f}s")
        print()
        rows.append({
            'metodoa': method,
            'supergelaxka_tamaina': size,
            'denbora_s': round(total, 2),
        })

os.makedirs(os.path.dirname(OUTPUT), exist_ok=True)
with open(OUTPUT, 'w', newline='') as f:
    writer = csv.DictWriter(f, fieldnames=['metodoa', 'supergelaxka_tamaina', 'denbora_s'])
    writer.writeheader()
    writer.writerows(rows)

print(f"\nSaved: {OUTPUT}")
print(f"\n{'Method':<8} {'Size':<8} {'Time (s)':>12}")
print("-" * 32)
for r in rows:
    print(f"{r['metodoa']:<8} {r['supergelaxka_tamaina']:<8} {r['denbora_s']:>12.1f}")
