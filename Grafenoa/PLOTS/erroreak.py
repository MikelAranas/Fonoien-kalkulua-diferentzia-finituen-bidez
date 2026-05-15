import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

def load_dat(filename):
    data = np.loadtxt(filename)
    s = data[:, 0]
    freqs = np.abs(data[:, 1:])
    return s, freqs

sizes   = ['2x2x1', '4x4x1', '7x7x1']
methods = ['DFPT','DF']
ref_size = '7x7x1'

data = {}
for method in methods:
    data[method] = {}
    for size in sizes:
        fname = f'{method}_{size}.dat'
        try:
            data[method][size] = load_dat(fname)
            print(f"Kargatuta: {fname}")
        except FileNotFoundError:
            print(f"Ez aurkitu: {fname}")
            data[method][size] = None

def relative_error(s_ref, f_ref, s_i, f_i):
    umbral = 10.0  # cm-1, ignorar modos muy pequeños (ZA cerca de Gamma)
    errors = []
    for b in range(f_ref.shape[1]):
        f_i_interp = np.interp(s_ref, s_i, f_i[:, b])
        mask = np.abs(f_ref[:, b]) > umbral
        if mask.sum() == 0:
            continue
        err = np.abs(f_ref[mask, b] - f_i_interp[mask]) / np.abs(f_ref[mask, b])
        errors.append(err)
    return np.concatenate(errors)

print("\n" + "="*60)
print(f"{'Metodo':<8} {'Tamaina':<8} {'Errore batez best. (%)':>22} {'Errore max. (%)':>16}")
print("="*60)

results = {method: {'sizes': [], 'mean': [], 'max': []} for method in methods}
size_to_N = {'2x2x1': 2, '3x3x1': 3, '4x4x1': 4, '6x6x1': 6}

for method in methods:
    ref = data[method][ref_size]
    if ref is None:
        continue
    s_ref, f_ref = ref
    for size in sizes[:-1]:
        if data[method][size] is None:
            continue
        s_i, f_i = data[method][size]
        errors = relative_error(s_ref, f_ref, s_i, f_i)
        mean_err = errors.mean() * 100
        max_err  = errors.max()  * 100
        print(f"{method:<8} {size:<8} {mean_err:>22.3f} {max_err:>16.3f}")
        results[method]['sizes'].append(size)
        results[method]['mean'].append(mean_err)
        results[method]['max'].append(max_err)

colors = {'DFPT': '#1f77b4', 'DF': '#d62728'}

fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('Konbergentzia: errore erlatiboa 6×6×1 erreferentziarekiko — Grafenoa', fontsize=13)

for ax, key, ylabel in zip(axes,
                            ['mean', 'max'],
                            ['Errore erlatiboa batez bestekoa (%)', 'Errore erlatiboa maximoa (%)']):
    for method in methods:
        if not results[method]['sizes']:
            continue
        Ns = [size_to_N[s] for s in results[method]['sizes']]
        vals = results[method][key]
        ax.plot(Ns, vals, 'o-', color=colors[method], linewidth=2,
                markersize=8, label=method)
    ax.set_xlabel('Supergelaxkaren tamaina N', fontsize=12)
    ax.set_ylabel(ylabel, fontsize=12)
    ax.set_xticks([2, 3, 4])
    ax.set_xticklabels(['2×2×1', '3×3×1', '4×4×1'])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('errore_erlatiboa_grafeno.pdf', bbox_inches='tight')
print("\nGordeta: errore_erlatiboa_grafeno.pdf")
