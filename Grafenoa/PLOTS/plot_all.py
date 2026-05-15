import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def load_dat(filename):
    data = np.loadtxt(filename)
    s = data[:, 0]
    freqs = np.abs(data[:, 1:])  # forzar positivo
    return s, freqs

sizes = ['2x2x1', '4x4x1', '7x7x1']
sym_s      = [0.000000, 0.577350, 0.910683, 1.577349]
sym_labels = [r'$\Gamma$', 'M', 'K', r'$\Gamma$']

# ============================================================
# FIGURA 1: Convergencia DFPT
# ============================================================
size_colors = {'2x2x1': '#1f77b4', '4x4x1': '#2ca02c', '7x7x1': 'black'}
size_lw     = {'2x2x1': 1.0, '4x4x1': 1.0, '7x7x1': 2.0}

fig, ax = plt.subplots(figsize=(8, 6))

for size in sizes:
    fname = f'DFPT_{size}.dat'
    try:
        s, freqs = load_dat(fname)
    except FileNotFoundError:
        print(f"No encontrado: {fname}")
        continue
    color = size_colors[size]
    lw = size_lw[size]
    for i in range(freqs.shape[1]):
        ax.plot(s, freqs[:, i], color=color, linewidth=lw,
                label=size if i == 0 else '_nolegend_')

for ss in sym_s:
    ax.axvline(x=ss, color='gray', linewidth=0.7, linestyle='--', alpha=0.6)
for ss, label in zip(sym_s, sym_labels):
    ax.text(ss, -25, label, ha='center', va='top', fontsize=13)

ax.set_xlim(sym_s[0], sym_s[-1])
ax.set_ylim(-5, 1700)
ax.set_ylabel(r'$\omega$ (cm$^{-1}$)', fontsize=13)
ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
ax.tick_params(axis='y', labelsize=12)
ax.axhline(y=0, color='black', linewidth=0.5)

legend_handles = [
    mlines.Line2D([], [], color=size_colors['2x2x1'], linewidth=1.5, label='2×2×1'),
    mlines.Line2D([], [], color=size_colors['4x4x1'], linewidth=1.5, label='4×4×1'),
    mlines.Line2D([], [], color=size_colors['7x7x1'], linewidth=2.5, label='7×7×1 (erreferentzia)'),
]
ax.legend(handles=legend_handles, fontsize=11)
ax.set_title('Grafenoa — DFPT konbergentzia', fontsize=13)
plt.tight_layout()
plt.savefig('konbergentzia_grafeno.pdf', bbox_inches='tight')
print("Gordeta: konbergentzia_grafeno.pdf")
plt.close()

# ============================================================
# FIGURA 2: DF (puntos) vs DFPT (líneas) per size
# ============================================================
fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for ax, size in zip(axes, sizes):
    try:
        s_dfpt, f_dfpt = load_dat(f'DFPT_{size}.dat')
    except FileNotFoundError:
        ax.set_title(f'{size} — DFPT ez aurkitu')
        continue
    try:
        s_df, f_df = load_dat(f'DF_{size}.dat')
        has_df = True
    except FileNotFoundError:
        has_df = False

    for i in range(f_dfpt.shape[1]):
        ax.plot(s_dfpt, f_dfpt[:, i], color='#1f77b4', linewidth=1.5,
                label='DFPT' if i == 0 else '_nolegend_')
    if has_df:
        step = max(1, len(s_df) // 40)
        for i in range(f_df.shape[1]):
            ax.plot(s_df[::step], f_df[::step, i], 'o', color='#d62728',
                    markersize=3.5, markeredgewidth=0.3, markeredgecolor='white',
                    label='DF' if i == 0 else '_nolegend_')

    for ss in sym_s:
        ax.axvline(x=ss, color='gray', linewidth=0.7, linestyle='--', alpha=0.6)
    for ss, label in zip(sym_s, sym_labels):
        ax.text(ss, -25, label, ha='center', va='top', fontsize=11)

    all_freqs = f_dfpt.flatten()
    if has_df:
        all_freqs = np.concatenate([all_freqs, f_df.flatten()])
    ax.set_xlim(sym_s[0], sym_s[-1])
    ax.set_ylim(-5, all_freqs.max() * 1.05)
    ax.set_title(size, fontsize=13)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=11)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel(r'$\omega$ (cm$^{-1}$)', fontsize=11)
    ax.legend(fontsize=10)

fig.suptitle('DF (puntuak) vs DFPT (lerroak) — Grafenoa', fontsize=14)
plt.tight_layout()
plt.savefig('DF_vs_DFPT_grafeno.pdf', bbox_inches='tight')
print("Gordeta: DF_vs_DFPT_grafeno.pdf")
