#!/usr/bin/env python3
import os
import sys

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "modules"))
import marrazketa as mz

cfg = {
    "dir": HERE,
    "sizes": ["2x2x2", "3x3x3", "4x4x4", "6x6x6"],
    "sym_s": [0.000, 0.707107, 0.957107, 1.707107, 2.319479],
    "sym_labels": [r"$\Gamma$", r"$X$", r"$K$", r"$\Gamma$", r"$L$"],
    "unit": r"$\omega$ (cm$^{-1}$)",
    "force_abs": False,
    "panel_grid": (2, 2),
    "figsize": (10.0, 7.2),
    "N_values": [2, 3, 4, 6],
}
NAME = "al"
CSV = os.path.join(HERE, "denborak.csv")
Q_IRR = {"2x2x2": 3, "3x3x3": 4, "4x4x4": 8, "6x6x6": 16}
N_SCF = 6

if __name__ == "__main__":
    data = mz.load_all(cfg)
    mz.fig_methods(HERE, NAME, cfg, data)
    mz.fig_convergence(HERE, NAME, cfg, data)
    mz.fig_errors(HERE, NAME, cfg, data)
    mz.fig_denbora(HERE, NAME, CSV)
    mz.emaitzak_taula(HERE, NAME, cfg, data, CSV, N_SCF, Q_IRR)
