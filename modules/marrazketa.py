#!/usr/bin/env python3

import csv
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

plt.rcParams.update({
    "font.family": "serif",
    "mathtext.fontset": "cm",
    "font.size": 11,
    "axes.titlesize": 12,
    "axes.labelsize": 12,
    "legend.fontsize": 10,
    "xtick.labelsize": 11,
    "ytick.labelsize": 10,
    "axes.linewidth": 0.8,
    "figure.constrained_layout.use": True,
})

C_DFPT = "#1f77b4"
C_DF = "#d62728"
C_HSE = "#ff7f0e"
C_EXP = "black"
SIZE_COLORS = ["#1f77b4", "#ff7f0e", "#2ca02c", "black"]

def load_dat(path, force_abs=False):
    if not os.path.isfile(path):
        return None
    data = np.loadtxt(path)
    s, freqs = data[:, 0], data[:, 1:]
    if force_abs:
        freqs = np.abs(freqs)
    return s, freqs

def load_all(cfg):
    out = {}
    for method in ("DFPT", "DF"):
        out[method] = {}
        for size in cfg["sizes"]:
            path = os.path.join(cfg["dir"], f"{method}_{size}.dat")
            out[method][size] = load_dat(path, cfg["force_abs"])
    return out

def load_exp(cfg, s_max):
    path = os.path.join(cfg["dir"], "EXP.csv")
    data = np.loadtxt(path, delimiter=",")
    x, y = data[:, 0], data[:, 1]
    return x * (s_max / x.max()), y

def si_sym_points(azterketa_dir):
    b1 = np.array([-0.707107, -0.707107, 0.707107])
    b2 = np.array([0.707107, 0.707107, 0.707107])
    b3 = np.array([-0.707107, 0.707107, -0.707107])
    B = np.vstack([b1, b2, b3])
    sym_pts = [
        (r"$\Gamma$", np.array([0.0, 0.0, 0.0])),
        (r"$X$", np.array([0.50, 0.50, 0.00])),
        (r"$K$", np.array([0.375, 0.625, 0.00])),
        (r"$\Gamma$", np.array([0.0, 1.0, 0.0])),
        (r"$L$", np.array([0.0, 1.5, 0.0])),
    ]
    q_path = np.loadtxt(os.path.join(azterketa_dir, "path.dat"), usecols=(0, 1, 2))
    q_cart = q_path @ B
    s_path = np.zeros(len(q_path))
    for i in range(1, len(q_path)):
        s_path[i] = s_path[i - 1] + np.linalg.norm(q_cart[i] - q_cart[i - 1])

    def s_at_q(q_target):
        best_s, best_d2 = None, np.inf
        for i in range(len(q_path) - 1):
            a, b = q_cart[i], q_cart[i + 1]
            v = b - a
            vv = np.dot(v, v)
            if vv == 0:
                continue
            t = np.clip(np.dot(q_target @ B - a, v) / vv, 0, 1)
            d2 = np.sum((q_target @ B - (a + t * v)) ** 2)
            if d2 < best_d2:
                best_d2, best_s = d2, s_path[i] + t * np.linalg.norm(v)
        return best_s

    return [s_at_q(q) for _, q in sym_pts], [lab for lab, _ in sym_pts]

def s_at_q_si(azterketa_dir, q_target):
    b1 = np.array([-0.707107, -0.707107, 0.707107])
    b2 = np.array([0.707107, 0.707107, 0.707107])
    b3 = np.array([-0.707107, 0.707107, -0.707107])
    B = np.vstack([b1, b2, b3])
    q_path = np.loadtxt(os.path.join(azterketa_dir, "path.dat"), usecols=(0, 1, 2))
    q_cart = q_path @ B
    s_path = np.zeros(len(q_path))
    for i in range(1, len(q_path)):
        s_path[i] = s_path[i - 1] + np.linalg.norm(q_cart[i] - q_cart[i - 1])
    qc = q_target @ B
    best_s, best_d2 = None, np.inf
    for i in range(len(q_path) - 1):
        a, b = q_cart[i], q_cart[i + 1]
        v = b - a
        vv = np.dot(v, v)
        if vv == 0:
            continue
        t = np.clip(np.dot(qc - a, v) / vv, 0, 1)
        d2 = np.sum((qc - (a + t * v)) ** 2)
        if d2 < best_d2:
            best_d2, best_s = d2, s_path[i] + t * np.linalg.norm(v)
    return best_s

def deco_axis(ax, cfg, ymax, ymin=0.0):
    for ss in cfg["sym_s"]:
        ax.axvline(x=ss, color="gray", linewidth=0.6, linestyle="--",
                   alpha=0.6, zorder=0)
    ax.set_xticks(cfg["sym_s"])
    ax.set_xticklabels(cfg["sym_labels"])
    ax.tick_params(axis="x", length=0)
    ax.set_xlim(cfg["sym_s"][0], cfg["sym_s"][-1])
    ax.set_ylim(ymin, ymax)
    ax.set_ylabel(cfg["unit"])

def fig_methods(outdir, name, cfg, data, exp=None):
    nrow, ncol = cfg["panel_grid"]
    fig, axes = plt.subplots(nrow, ncol, figsize=cfg["figsize"], sharey=True)
    axes = np.atleast_1d(axes).flatten()

    ymax = 0.0
    for size in cfg["sizes"]:
        for method in ("DFPT", "DF"):
            if data[method][size] is not None:
                ymax = max(ymax, data[method][size][1].max())
    if exp is not None:
        ymax = max(ymax, exp[1].max())
    ymax *= 1.05

    for ax, size in zip(axes, cfg["sizes"]):
        s_dfpt, f_dfpt = data["DFPT"][size]
        for i in range(f_dfpt.shape[1]):
            ax.plot(s_dfpt, f_dfpt[:, i], color=C_DFPT, lw=1.4,
                    label="DFPT" if i == 0 else "_nolegend_")
        if data["DF"][size] is not None:
            s_df, f_df = data["DF"][size]
            step = max(1, len(s_df) // 45)
            for i in range(f_df.shape[1]):
                ax.plot(s_df[::step], f_df[::step, i], "o", color=C_DF,
                        ms=3.0, mew=0.3, mec="white",
                        label="DF" if i == 0 else "_nolegend_")
        if exp is not None:
            ax.scatter(exp[0], exp[1], s=16, marker="D", facecolors="none",
                       edgecolors=C_EXP, linewidths=0.7, zorder=5, label="Esp.")
        deco_axis(ax, cfg, ymax)
        if ax is not axes[0]:
            ax.set_ylabel("")
        ax.set_title(size.replace("x", r"$\times$"))
        ax.legend(loc="lower right" if name == "al" else "best", framealpha=0.9)

    for ax in axes[len(cfg["sizes"]):]:
        ax.set_visible(False)

    _save(fig, outdir, f"{name}_metodoak.pdf")

def fig_convergence(outdir, name, cfg, data):
    fig, ax = plt.subplots(figsize=(7.0, 5.0))
    sizes = cfg["sizes"]
    ymax = 0.0
    for k, size in enumerate(sizes):
        s, f = data["DFPT"][size]
        ymax = max(ymax, f.max())
        is_ref = size == sizes[-1]
        color = "black" if is_ref else SIZE_COLORS[k]
        lw = 2.0 if is_ref else 1.0
        lab = size.replace("x", r"$\times$") + (" (erref.)" if is_ref else "")
        for i in range(f.shape[1]):
            ax.plot(s, f[:, i], color=color, lw=lw,
                    label=lab if i == 0 else "_nolegend_",
                    zorder=3 if is_ref else 2)
    deco_axis(ax, cfg, ymax * 1.05)
    ax.legend()
    _save(fig, outdir, f"{name}_konbergentzia.pdf")

def deviations(data, method, sizes):
    ref = data[method][sizes[-1]]
    if ref is None:
        ref = data["DFPT"][sizes[-1]]
    s_ref, f_ref = ref
    means, maxima = [], []
    for size in sizes[:-1]:
        s_i, f_i = data[method][size]
        diffs = [np.abs(f_ref[:, b] - np.interp(s_ref, s_i, f_i[:, b]))
                 for b in range(f_ref.shape[1])]
        diffs = np.concatenate(diffs)
        means.append(diffs.mean())
        maxima.append(diffs.max())
    return np.array(means), np.array(maxima)

def fig_errors(outdir, name, cfg, data):
    sizes = cfg["sizes"]
    x = cfg["N_values"][:-1]
    unit = cfg["unit"].split("(")[1].rstrip(")")
    ref_lab = sizes[-1].replace("x", r"$\times$")
    dev = {m: deviations(data, m, sizes) for m in ("DFPT", "DF")}
    for eskala in ("log", "lineal"):
        fig, axes = plt.subplots(1, 2, figsize=(10.0, 4.0))
        for method, color, marker in (("DFPT", C_DFPT, "s"), ("DF", C_DF, "o")):
            mean, mx = dev[method]
            axes[0].plot(x, mean, marker=marker, color=color, label=method)
            axes[1].plot(x, mx, marker=marker, color=color, label=method)
        axes[0].set_ylabel(rf"$\langle|\Delta\omega|\rangle$ ({unit})")
        axes[1].set_ylabel(rf"$\max|\Delta\omega|$ ({unit})")
        for ax in axes:
            ax.set_xlabel(r"Supergelaxkaren tamaina $N$")
            ax.set_xticks(cfg["N_values"][:-1])
            ax.set_xticklabels([s.replace("x", r"$\times$") for s in sizes[:-1]])
            if eskala == "log":
                ax.set_yscale("log")
            ax.grid(alpha=0.3)
            ax.legend(title=f"Erref.: {ref_lab}")
        _save(fig, outdir, f"{name}_errorea_{eskala}.pdf")

def fig_za(outdir, name, cfg, data):
    fig, ax = plt.subplots(figsize=(7.0, 4.6))
    for k, size in enumerate(cfg["sizes"]):
        s, f = data["DFPT"][size]
        za = f.min(axis=1)
        is_ref = size == cfg["sizes"][-1]
        color = "black" if is_ref else SIZE_COLORS[k]
        lab = size.replace("x", r"$\times$") + (" (erref.)" if is_ref else "")
        ax.plot(s, za, color=color, lw=2.0 if is_ref else 1.2, label=lab,
                zorder=3 if is_ref else 2)
        if data["DF"][size] is not None:
            s_df, f_df = data["DF"][size]
            step = max(1, len(s_df) // 60)
            ax.plot(s_df[::step], f_df.min(axis=1)[::step], "o", color=color,
                    ms=2.8, mew=0.3, mec="white")
    deco_axis(ax, cfg, 700.0)
    ax.legend(title="Lerroak DFPT, puntuak DF")
    _save(fig, outdir, f"{name}_ZA.pdf")

def fig_si_hse(outdir, name, cfg, exp):
    s_pbe, f_pbe = load_dat(os.path.join(cfg["dir"], "DF_2x2x2.dat"))
    s_hse, f_hse = load_dat(os.path.join(cfg["dir"], "DF_2x2x2_hse.dat"))
    fig, ax = plt.subplots(figsize=(7.6, 5.2))
    for i in range(f_pbe.shape[1]):
        ax.plot(s_pbe, f_pbe[:, i], color=C_DFPT, lw=1.4,
                label="DF (PBE)" if i == 0 else "_nolegend_")
    for i in range(f_hse.shape[1]):
        ax.plot(s_hse, f_hse[:, i], color=C_HSE, lw=1.4, ls="--",
                label="DF (HSE)" if i == 0 else "_nolegend_")
    ax.scatter(exp[0], exp[1], s=18, marker="D", facecolors="none",
               edgecolors=C_EXP, linewidths=0.7, zorder=5, label="Esp.")
    ymax = max(f_pbe.max(), f_hse.max(), exp[1].max()) * 1.05
    deco_axis(ax, cfg, ymax)
    ax.legend(loc="lower right")
    _save(fig, outdir, f"{name}_hse.pdf")

def fig_denbora(outdir, name, csv_path):
    rows = list(csv.DictReader(open(csv_path)))
    series = {}
    for method in ("DF", "DFPT"):
        pts = sorted((int(r["supergelaxka_tamaina"].split("x")[0]), float(r["denbora_s"]))
                     for r in rows if r["metodoa"] == method)
        if pts:
            series[method] = list(zip(*pts))
    for eskala in ("log", "lineal"):
        fig, ax = plt.subplots(figsize=(6.0, 4.6))
        for method, color, marker in (("DF", C_DF, "o"), ("DFPT", C_DFPT, "s")):
            if method in series:
                xs, ys = series[method]
                ax.plot(xs, ys, marker=marker, color=color, label=method)
        if eskala == "log":
            ax.set_yscale("log")
        ax.set_xlabel(r"Supergelaxkaren tamaina $N$")
        ax.set_ylabel(r"$t$ (s)")
        ax.grid(alpha=0.3, which="both")
        ax.legend()
        _save(fig, outdir, f"{name}_denbora_{eskala}.pdf")

def emaitzak_taula(outdir, name, cfg, data, csv_path, n_scf, q_irr, sym=None):
    lines = []
    w = lines.append
    w("=" * 60)
    w(f"  {name.upper()} - AZTERKETAREN ZENBAKIAK")
    w("=" * 60)

    w("\n1. KOSTU KONPUTAZIONALA")
    t = {(r["metodoa"], r["supergelaxka_tamaina"]): float(r["denbora_s"])
         for r in csv.DictReader(open(csv_path))}
    w(f"   (DF: {n_scf} SCF desplazatu)")
    w(f"   {'tamaina':<8} {'q_irr':>5} {'t_DFPT(s)':>10} {'t_DF(s)':>12} {'faktorea':>9}")
    for size, nq in q_irr.items():
        td = t.get(("DFPT", size))
        tf = t.get(("DF", size))
        if td is None:
            continue
        if tf is None:
            w(f"   {size:<8} {nq:>5} {td:>10.1f} {'exekutatzen':>12} {'---':>9}")
        else:
            w(f"   {size:<8} {nq:>5} {td:>10.1f} {tf:>12.1f} {tf/td:>8.1f}x")

    w("\n2. DF vs DFPT ADOSTASUNA TAMAINA BEREAN (batez besteko |dw|)")
    unit = cfg["unit"].split("(")[1].rstrip(")")
    for size in cfg["sizes"]:
        d1, d2 = data["DFPT"][size], data["DF"][size]
        if d1 is None or d2 is None:
            continue
        s1, f1 = d1
        s2, f2 = d2
        diffs = np.concatenate([np.abs(f1[:, b] - np.interp(s1, s2, f2[:, b]))
                                for b in range(f1.shape[1])])
        w(f"   {size}: <|dw|> = {diffs.mean():.4f}, max = {diffs.max():.4f} {unit}")

    if sym:
        w("\n3. MAIZTASUNAK PUNTU SIMETRIKOETAN (THz)")
        dsets = {
            "DFPT 6x6x6": data["DFPT"].get("6x6x6"),
            "DF   6x6x6": data["DF"].get("6x6x6"),
            "HSE  2x2x2": load_dat(os.path.join(cfg["dir"], "DF_2x2x2_hse.dat")),
            "PBE  2x2x2": data["DF"].get("2x2x2"),
        }
        for pt, qf in sym.items():
            sp = s_at_q_si(cfg["dir"], np.array(qf))
            w(f"   {pt} (s={sp:.4f}):")
            for nm, ds in dsets.items():
                if ds is None:
                    w(f"     {nm}: (exekutatzen)")
                    continue
                s, f = ds
                idx = np.argmin(np.abs(s - sp))
                freqs = sorted(f[idx])
                uniq = []
                for fr in freqs:
                    if uniq and abs(fr - uniq[-1][0]) < 0.05:
                        uniq[-1] = ((uniq[-1][0] * uniq[-1][1] + fr) / (uniq[-1][1] + 1),
                                    uniq[-1][1] + 1)
                    else:
                        uniq.append((fr, 1))
                w(f"     {nm}: " + ", ".join(f"{fr:.2f}(x{m})" for fr, m in uniq))

    out = os.path.join(outdir, "azterketa_emaitzak.txt")
    open(out, "w").write("\n".join(lines) + "\n")
    print("Saved:", out)

def _save(fig, outdir, fname):
    out = os.path.join(outdir, fname)
    fig.savefig(out)
    plt.close(fig)
    print("Saved:", out)
