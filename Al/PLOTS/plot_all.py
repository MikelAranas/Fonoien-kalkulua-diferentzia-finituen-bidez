import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def load_dat(filename):
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1:]

sizes = ['2x2x2', '3x3x3', '4x4x4', '6x6x6']
sym_s      = [0.000, 0.707107, 0.957107, 1.707107, 2.319479]
sym_labels = [r'$\Gamma$', 'X', 'K', r'$\Gamma$', 'L']

fig, axes = plt.subplots(2, 2, figsize=(14, 10))
axes = axes.flatten()

for ax, size in zip(axes, sizes):
    try:
        s_dfpt, f_dfpt = load_dat(f'DFPT_{size}.dat')
    except FileNotFoundError:
        ax.set_title(f'{size} — DFPT no encontrado')
        continue
    try:
        s_df, f_df = load_dat(f'DF_{size}.dat')
        has_df = True
    except FileNotFoundError:
        has_df = False

    # DFPT: líneas
    for i in range(f_dfpt.shape[1]):
        ax.plot(s_dfpt, f_dfpt[:, i], color='#1f77b4', linewidth=1.5,
                label='DFPT' if i == 0 else '_nolegend_')

    # DF: puntos
    if has_df:
        step = max(1, len(s_df) // 40)
        for i in range(f_df.shape[1]):
            ax.plot(s_df[::step], f_df[::step, i], 'o', color='#d62728',
                    markersize=3.5, markeredgewidth=0.3, markeredgecolor='white',
                    label='DF' if i == 0 else '_nolegend_')

    # Puntos de alta simetría
    for ss in sym_s:
        ax.axvline(x=ss, color='gray', linewidth=0.7, linestyle='--', alpha=0.6)
    for ss, label in zip(sym_s, sym_labels):
        ax.text(ss, -18, label, ha='center', va='top', fontsize=11)

    all_freqs = f_dfpt.flatten()
    if has_df:
        all_freqs = np.concatenate([all_freqs, f_df.flatten()])
    ymax = all_freqs.max() * 1.05

    ax.set_xlim(s_dfpt[0], s_dfpt[-1])
    ax.set_ylim(-5, ymax)
    ax.set_title(f'{size}', fontsize=13)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=11)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel(r'$\omega$ (cm$^{-1}$)', fontsize=11)
    ax.legend(fontsize=10)

fig.suptitle('DF (puntos) vs DFPT (líneas) — Al', fontsize=14)
plt.tight_layout()
plt.savefig('DF_vs_DFPT_per_size.pdf', bbox_inches='tight')
print("Gordeta: DF_vs_DFPT_per_size.pdf")
