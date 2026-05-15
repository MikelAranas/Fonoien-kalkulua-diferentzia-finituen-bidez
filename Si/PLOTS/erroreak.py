import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines

# --- Cargar datos ---
def load_dat(filename):
    data = np.loadtxt(filename)
    s = data[:, 0]
    freqs = data[:, 1:]
    return s, freqs

sizes = ['2x2x2', '3x3x3', '4x4x4', '6x6x6']
methods = ['DFPT','DF']

data = {}
for method in methods:
    data[method] = {}
    for size in sizes:
        fname = f'{method}_{size}.dat'
        try:
            data[method][size] = load_dat(fname)
            print(f"Cargado: {fname}")
        except FileNotFoundError:
            print(f"No encontrado: {fname}")
            data[method][size] = None

# Referencia: 6x6x6
ref_size = '6x6x6'

# --- Calcular error relativo ---
def relative_error(s_ref, f_ref, s_i, f_i):
    """
    Interpola f_i en los puntos s_ref y calcula el error relativo
    error(q) = |f_ref(q) - f_i(q)| / |f_ref(q)|
    Solo considera frecuencias > umbral para evitar division por cero cerca de Gamma
    """
    umbral = 5.0  # cm-1, ignorar frecuencias muy pequeñas
    errors = []
    n_branches = f_ref.shape[1]
    for b in range(n_branches):
        # Interpolar f_i en los puntos de s_ref
        f_i_interp = np.interp(s_ref, s_i, f_i[:, b])
        mask = np.abs(f_ref[:, b]) > umbral
        if mask.sum() == 0:
            continue
        err = np.abs(f_ref[mask, b] - f_i_interp[mask]) / np.abs(f_ref[mask, b])
        errors.append(err)
    errors = np.concatenate(errors)
    return errors

# --- Calcular y mostrar errores ---
print("\n" + "="*60)
print(f"{'Método':<8} {'Tamaño':<8} {'Error medio (%)':>16} {'Error max (%)':>14}")
print("="*60)

results = {method: {'sizes': [], 'mean': [], 'max': []} for method in methods}

for method in methods:
    ref = data[method][ref_size]
    if ref is None:
        print(f"{method}: referencia 6x6x6 no disponible")
        continue
    s_ref, f_ref = ref

    for size in sizes[:-1]:  # excluir la referencia
        if data[method][size] is None:
            continue
        s_i, f_i = data[method][size]
        errors = relative_error(s_ref, f_ref, s_i, f_i)
        mean_err = errors.mean() * 100
        max_err  = errors.max()  * 100
        print(f"{method:<8} {size:<8} {mean_err:>16.3f} {max_err:>14.3f}")
        results[method]['sizes'].append(size)
        results[method]['mean'].append(mean_err)
        results[method]['max'].append(max_err)

# --- Mapear tamaños a N ---
size_to_N = {'2x2x2': 2, '3x3x3': 3, '4x4x4': 4, '6x6x6': 6}
colors = {'DFPT': '#1f77b4', 'DF': '#d62728'}

# ============================================================
# FIGURA 1: Error medio y máximo vs N
# ============================================================
fig, axes = plt.subplots(1, 2, figsize=(13, 5))
fig.suptitle('Konbergentzia: errore erlatiboa 6×6×6 erreferentziarekiko', fontsize=13)

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
    ax.set_xticklabels(['2×2×2', '3×3×3', '4×4×4'])
    ax.legend(fontsize=11)
    ax.grid(True, alpha=0.3)

plt.tight_layout()
plt.savefig('errore_erlatiboa.pdf', bbox_inches='tight')
print("\nGordeta: errore_erlatiboa.pdf")

# ============================================================
# FIGURA 2: Dispersión — todas las curvas + referencia destacada
# ============================================================
# Puntos de alta simetria Si: W->X->G->L->K
# Calcular s_path en los puntos de simetria usando el DFPT 6x6x6 como referencia
b1 = np.array([-0.707107, -0.707107,  0.707107])
b2 = np.array([ 0.707107,  0.707107,  0.707107])
b3 = np.array([-0.707107,  0.707107, -0.707107])
B = np.vstack([b1, b2, b3])

sym_pts = [
    (r'$W$',      np.array([0.750, 0.500, 0.250])),
    (r'$X$',      np.array([0.500, 0.500, 0.000])),
    (r'$\Gamma$', np.array([0.000, 0.000, 0.000])),
    (r'$L$',      np.array([0.500, 0.500, 0.500])),
    (r'$K$',      np.array([0.750, 0.375, 0.375])),
]

def s_at_q(q_path, s_path, q_target):
    best_s, best_d2 = None, np.inf
    for i in range(len(q_path)-1):
        a, b = q_path[i], q_path[i+1]
        v = b - a
        vv = np.dot(v, v)
        if vv == 0: continue
        t = np.clip(np.dot(q_target - a, v) / vv, 0, 1)
        d2 = np.dot(q_target - (a + t*v), q_target - (a + t*v))
        if d2 < best_d2:
            best_d2 = d2
            best_s = s_path[i] + t * np.linalg.norm(v)
    return best_s

# Cargar path del DFPT 6x6x6 para calcular s_sym
q_path_ref = np.loadtxt('path.dat', usecols=(0,1,2))
q_cart_ref = q_path_ref @ B
s_ref_path = np.zeros(len(q_path_ref))
for i in range(1, len(q_path_ref)):
    s_ref_path[i] = s_ref_path[i-1] + np.linalg.norm(q_cart_ref[i] - q_cart_ref[i-1])

sym_s      = [s_at_q(q_path_ref, s_ref_path, qpt) for _, qpt in sym_pts]
sym_labels = [label for label, _ in sym_pts]
size_colors = {'2x2x2': '#1f77b4', '3x3x3': '#ff7f0e',
               '4x4x4': '#2ca02c', '6x6x6': 'black'}
size_lw     = {'2x2x2': 1.0, '3x3x3': 1.0, '4x4x4': 1.0, '6x6x6': 2.0}

fig, axes = plt.subplots(1, 2, figsize=(14, 6), sharey=True)
fig.subplots_adjust(wspace=0.05)

for ax, method in zip(axes, methods):
    if data[method][ref_size] is None:
        ax.set_title(f'{method} (sin datos)')
        continue

    for size in sizes:
        if data[method][size] is None:
            continue
        s, freqs = data[method][size]
        color = size_colors[size]
        lw = size_lw[size]
        for i in range(freqs.shape[1]):
            ax.plot(s, freqs[:, i], color=color, linewidth=lw,
                    label=size if i == 0 else '_nolegend_')

    for ss in sym_s:
        ax.axvline(x=ss, color='gray', linewidth=0.7, linestyle='--', alpha=0.6)
    for ss, label in zip(sym_s, sym_labels):
        ax.text(ss, -18, label, ha='center', va='top', fontsize=13)

    s_ref, f_ref = data[method][ref_size]
    ax.set_xlim(s_ref[0], s_ref[-1])
    # Calcular ylim dinámicamente
    all_freqs = np.concatenate([data[method][s][1].flatten() 
                                for s in sizes if data[method][s] is not None])
    ymax = all_freqs.max() * 1.05
    ax.set_ylim(-5, ymax)
    ax.set_title(method, fontsize=14)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=12)
    ax.axhline(y=0, color='black', linewidth=0.5)

axes[0].set_ylabel(r'Maiztasuna (cm$^{-1}$)', fontsize=13)

legend_handles = [
    mlines.Line2D([], [], color=size_colors['2x2x2'], linewidth=1.5, label='2×2×2'),
    mlines.Line2D([], [], color=size_colors['3x3x3'], linewidth=1.5, label='3×3×3'),
    mlines.Line2D([], [], color=size_colors['4x4x4'], linewidth=1.5, label='4×4×4'),
    mlines.Line2D([], [], color=size_colors['6x6x6'], linewidth=2.5, label='6×6×6 (erreferentzia)'),
]
fig.legend(handles=legend_handles, loc='upper center', ncol=4,
           fontsize=11, frameon=True, bbox_to_anchor=(0.5, 1.02))

plt.savefig('konbergentzia_6x6x6.pdf', bbox_inches='tight')
print("Gordeta: konbergentzia_6x6x6.pdf")
