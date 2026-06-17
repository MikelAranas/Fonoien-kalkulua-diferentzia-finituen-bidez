#!/usr/bin/env python3
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "modules"))
import marrazketa as mz

cfg = {
    "dir": HERE,
    "sizes": ["2x2x1", "4x4x1", "7x7x1"],
    "sym_s": [0.000000, 0.577350, 0.910683, 1.577349],
    "sym_labels": [r"$\Gamma$", r"$M$", r"$K$", r"$\Gamma$"],
    "unit": r"$\omega$ (cm$^{-1}$)",
    "force_abs": True,
    "panel_grid": (1, 3),
    "figsize": (12.0, 4.0),
    "N_values": [2, 4, 7],
}
NAME = "grafenoa"
CSV = os.path.join(HERE, "denborak.csv")
Q_IRR = {"2x2x1": 2, "4x4x1": 4, "7x7x1": 8}
N_SCF = 12

if __name__ == "__main__":
    data = mz.load_all(cfg)
    mz.fig_methods(HERE, NAME, cfg, data)
    mz.fig_convergence(HERE, NAME, cfg, data)
    mz.fig_errors(HERE, NAME, cfg, data)
    mz.fig_za(HERE, NAME, cfg, data)
    mz.fig_denbora(HERE, NAME, CSV)
    mz.emaitzak_taula(HERE, NAME, cfg, data, CSV, N_SCF, Q_IRR)
