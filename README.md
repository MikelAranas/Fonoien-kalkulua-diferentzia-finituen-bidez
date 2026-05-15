# Fonoien kalkulua diferentzia finituen bidez

Gradu Amaierako Lana — Fisikako Gradua (UPV/EHU)

Biltegi honetan hiru materialen fonoi-bandak kalkulatzen dira
(Si, Al, grafenoa), bi metodo erabiliz:

- **DFPT** — *Density Functional Perturbation Theory*, Quantum ESPRESSO bidez.
- **Diferentzia finituak (inplementazio propioa)** — Python kode propioa
  (`numpy` + `matplotlib`) + QE.

## Egitura

```
TFG/
├── modules/
│   └── svecs_module.py
├── pseudopotentzialak/
│   └── {materiala}.upf
└── {materiala}/                       materiala ∈ { Si, Al, Grafenoa }
    ├── DF/{N1xN2xN3}/
    ├── DFPT/{N1xN2xN3}/
    ├── PLOTS/
    ├── denborak.csv
    ├── denborak.py
    ├── grafikoa_denb.py
    └── escalado_tiempos.pdf
```

