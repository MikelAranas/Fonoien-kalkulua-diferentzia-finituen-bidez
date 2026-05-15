# Aluminio (Al)

Metal monoatómico, estructura FCC. Solo tiene **modos acústicos** (3) por
celda primitiva — todos cero en Γ por invarianza traslacional.

## Parámetros

| Magnitud | Valor |
|----------|-------|
| Constante de red `a` | 4.05 Å |
| Base | Al en (0, 0, 0) (1 átomo por celda primitiva) |
| Masa atómica | 26.9815 uma |
| Pseudopotencial | `Al.upf` — ATOMPAW, **PBE** |
| Cutoff (ecutwfc) | 38 Ry |
| Ocupaciones | `smearing = 'mp'` (Methfessel–Paxton), `degauss = 0.02` Ry |
| `conv_thr` (SCF) | 10⁻¹⁶ Ry |
| DELTA (desplazamiento DF) | 0.01 Å |
| Tamaños de supercelda DF | 2×2×2, 3×3×3, 4×4×4, 6×6×6 |
| k-grid de la celda primitiva | 12×12×12 |
| Malla q DFPT | 6×6×6 |
| Camino q | Γ → X → K → Γ → L |

Al ser metal, el k-grid es relativamente denso (12×12×12 en la celda
primitiva) y se aplica smearing de Methfessel–Paxton para que el
muestreo de la superficie de Fermi converja correctamente.

## Resultados — frecuencia máxima (LA, frontera de zona, cm⁻¹)

| N | DF | DFPT | Δ (cm⁻¹) |
|---|-----|------|-----------|
| 2 | 318.51 | 318.55 | −0.04 |
| 3 | 324.05 | 324.01 | +0.04 |
| 4 | 318.51 | 318.55 | −0.04 |
| 6 | **318.50** | **318.56** | **−0.06** |

DF y DFPT coinciden al nivel de la centésima de cm⁻¹. La comparación
completa de bandas está en `PLOTS/DF_vs_DFPT_per_size.pdf`.

## Notas específicas

- **DF y DFPT usan el mismo pseudopotencial PBE** (un único archivo en
  `pseudopotentzialak/Al.upf`).
- Solo hay **6 SCFs desplazados** por supercelda (1 átomo de la base × 6
  direcciones), frente a los 12 de Si y grafeno.
- El orden de los átomos en `sc_positions` (en `D.py`) sigue el mismo
  patrón que el bloque `ATOMIC_POSITIONS` del input QE: `ix` exterior,
  `iz` interior. Durante la unificación del módulo `svecs_module.py`
  detectamos que una versión anterior usaba un orden distinto generado
  automáticamente; las bandas salían correctas por la alta simetría
  FCC, pero el etiquetado era inconsistente. Corregido en `D.py` actual.

## Reproducir

```bash
cd DF/6x6x6/out
python IFC_murriztue.py
python D.py
bash konparaketa.sh

cd ../../../PLOTS
python datuak_batu.py
python plot_all.py
python erroreak.py

cd ..
python denborak.py
python grafikoa_denb.py
```
