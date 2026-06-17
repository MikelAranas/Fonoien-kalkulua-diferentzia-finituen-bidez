#!/usr/bin/env python3
import os
import sys

import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__))
sys.path.insert(0, os.path.join(HERE, "..", "..", "modules"))
import marrazketa as mz

sym_s, sym_labels = mz.si_sym_points(HERE)
cfg = {
    "dir": HERE,
    "sizes": ["2x2x2", "3x3x3", "4x4x4", "6x6x6"],
    "sym_s": sym_s,
    "sym_labels": sym_labels,
    "unit": r"$\omega$ (THz)",
    "force_abs": False,
    "panel_grid": (2, 2),
    "figsize": (10.0, 7.2),
    "N_values": [2, 3, 4, 6],
}
NAME = "si"
CSV = os.path.join(HERE, "denborak.csv")
Q_IRR = {"2x2x2": 3, "3x3x3": 4, "4x4x4": 8, "6x6x6": 16}
N_SCF = 12
SYM = {
    "Gamma": [0.0, 1.0, 0.0],
    "X": [0.50, 0.50, 0.0],
    "L": [0.0, 1.5, 0.0],
}

if __name__ == "__main__":
    data = mz.load_all(cfg)
    s_ref = data["DFPT"]["6x6x6"][0]
    exp = mz.load_exp(cfg, s_ref[-1])
    mz.fig_methods(HERE, NAME, cfg, data, exp=exp)
    mz.fig_convergence(HERE, NAME, cfg, data)
    mz.fig_errors(HERE, NAME, cfg, data)
    mz.fig_si_hse(HERE, NAME, cfg, exp)
    mz.fig_denbora(HERE, NAME, CSV)
    mz.emaitzak_taula(HERE, NAME, cfg, data, CSV, N_SCF, Q_IRR, sym=SYM)
