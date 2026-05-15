"""
Read denborak.csv and plot DF vs DFPT timing vs supercell size.
"""
import csv
import re

import matplotlib.pyplot as plt


CSV_FILE   = 'denborak.csv'
OUTPUT_PDF = 'escalado_tiempos.pdf'


def parse_size(label):
    """e.g. '2x2x2' → 2,  '7x7x1' → 7."""
    return int(re.match(r'(\d+)', label).group(1))


with open(CSV_FILE) as f:
    rows = list(csv.DictReader(f))

times = {'DF': [], 'DFPT': []}
sizes = {'DF': [], 'DFPT': []}
labels = {'DF': [], 'DFPT': []}
for row in rows:
    method = row['metodoa']
    size   = parse_size(row['supergelaxka_tamaina'])
    t      = float(row['denbora_s'])
    sizes[method].append(size)
    times[method].append(t)
    labels[method].append(row['supergelaxka_tamaina'])

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(sizes['DF'],   times['DF'],   'o-',
        color='#d62728', lw=2, ms=8, label='DF')
ax.plot(sizes['DFPT'], times['DFPT'], 's-',
        color='#1f77b4', lw=2, ms=8, label='DFPT')
ax.set_xlabel('Supergelaxken tamaina', fontsize=13)
ax.set_ylabel('t (s)', fontsize=13)
ax.set_title('DFPT vs DF', fontsize=14)
ax.set_xticks(sizes['DF'])
ax.set_xticklabels(labels['DF'])
ax.legend(fontsize=12)
ax.grid(True, alpha=0.3)
plt.tight_layout()
plt.savefig(OUTPUT_PDF, bbox_inches='tight')
print(f"Guardado: {OUTPUT_PDF}")
