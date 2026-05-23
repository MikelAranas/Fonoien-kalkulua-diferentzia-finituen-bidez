import numpy as np
import matplotlib.pyplot as plt
import matplotlib.lines as mlines
import os

def load_dat(filename):
    data = np.loadtxt(filename)
    return data[:, 0], data[:, 1:]

def load_exp(filename):
    """Carga EXP.csv: columna 0 = x en el path, columna 1 = frecuencia (THz)"""
    data = np.loadtxt(filename, delimiter=',')
    return data[:, 0], data[:, 1]

sizes = ['2x2x2', '3x3x3', '4x4x4', '6x6x6']

b1 = np.array([-0.707107, -0.707107,  0.707107])
b2 = np.array([ 0.707107,  0.707107,  0.707107])
b3 = np.array([-0.707107,  0.707107, -0.707107])
B = np.vstack([b1, b2, b3])

sym_pts = [
    (r'$\Gamma$', np.array([0.0,  0.0,  0.0])),
    (r'$X$',      np.array([0.50, 0.50, 0.00])),
    (r'$K$',      np.array([0.375, 0.625, 0.00])),
    (r'$\Gamma*$',np.array([0.0,  1.0,  0.0])),
    (r'$L$',      np.array([0.0,  1.5,  0.0])),
]

def s_at_q(q_path, s_path, q_target):
    best_s, best_d2 = None, np.inf
    for i in range(len(q_path) - 1):
        a, b = q_path[i], q_path[i + 1]
        v = b - a
        vv = np.dot(v, v)
        if vv == 0:
            continue
        t = np.clip(np.dot(q_target - a, v) / vv, 0, 1)
        d2 = np.dot(q_target - (a + t * v), q_target - (a + t * v))
        if d2 < best_d2:
            best_d2 = d2
            best_s = s_path[i] + t * np.linalg.norm(v)
    return best_s

# Cargar path del DFPT 6x6x6
q_path_ref = np.loadtxt('path.dat', usecols=(0, 1, 2))
q_cart_ref = q_path_ref @ B
s_ref_path = np.zeros(len(q_path_ref))
for i in range(1, len(q_path_ref)):
    s_ref_path[i] = s_ref_path[i - 1] + np.linalg.norm(q_cart_ref[i] - q_cart_ref[i - 1])

sym_s      = [s_at_q(q_path_ref, s_ref_path, qpt) for _, qpt in sym_pts]
sym_labels = [label for label, _ in sym_pts]

# Cargar datos experimentales (EXP.csv en el mismo directorio)
exp_file = 'EXP.csv'
has_exp = os.path.isfile(exp_file)
if has_exp:
    x_exp_raw, y_exp = load_exp(exp_file)
    # Escalar x_exp al rango del path de referencia
    # Los datos experimentales tienen x en [0, ~2.32], igual que s_ref_path
    # Si el rango del path de referencia es diferente, escalar proporcionalmente
    x_max_ref = s_ref_path[-1]
    x_max_exp = x_exp_raw.max()
    x_exp = x_exp_raw * (x_max_ref / x_max_exp)
else:
    print(f"[AVISO] No se encontró {exp_file} — se plotea sin datos experimentales")

fig, axes = plt.subplots(2, 2, figsize=(18, 14), constrained_layout=True)
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

    # Experimental: puntos grises con borde
    if has_exp:
        ax.scatter(x_exp, y_exp,
                   s=18, color='#2ca02c', marker='s',
                   linewidths=0.3, edgecolors='white',
                   zorder=5, label='Exp.')

    # Puntos de alta simetría
    for ss in sym_s:
        ax.axvline(x=ss, color='gray', linewidth=0.7, linestyle='--', alpha=0.6)
    for ss, label in zip(sym_s, sym_labels):
        ax.text(ss, -18, label, ha='center', va='top', fontsize=11)

    all_freqs = f_dfpt.flatten()
    if has_df:
        all_freqs = np.concatenate([all_freqs, f_df.flatten()])
    if has_exp:
        all_freqs = np.concatenate([all_freqs, y_exp])
    ymax = all_freqs.max() * 1.05

    ax.set_xlim(s_dfpt[0], s_dfpt[-1])
    ax.set_ylim(0.0, ymax)
    ax.set_title(f'{size}', fontsize=13)
    ax.tick_params(axis='x', which='both', bottom=False, labelbottom=False)
    ax.tick_params(axis='y', labelsize=11)
    ax.axhline(y=0, color='black', linewidth=0.5)
    ax.set_ylabel(r'$\omega$ (THz)', fontsize=11)
    ax.legend(fontsize=10)

fig.suptitle('DF (puntos) vs DFPT (líneas) — Si', fontsize=14)
plt.savefig('DF_vs_DFPT_per_size.pdf', bbox_inches='tight')
print("Gordeta: DF_vs_DFPT_per_size.pdf")