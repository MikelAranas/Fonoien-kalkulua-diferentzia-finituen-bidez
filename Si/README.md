# Silicio (Si)

Aislante covalente, estructura diamante (FCC con base de 2 átomos).
Comparación DF / DFPT en cuatro tamaños de supercelda.

## Parámetros

| Magnitud | Valor |
|----------|-------|
| Constante de red `a` | 5.43 Å |
| Base | Si en (0, 0, 0), Si en (¼, ¼, ¼) |
| Masa atómica | 28.086 uma |
| Pseudopotencial | `Si.upf` — ONCVPSP, **LDA** |
| Cutoff (ecutwfc) | 35 Ry |
| Ocupaciones | `fixed` (aislante) |
| `conv_thr` (SCF) | 10⁻¹² Ry |
| DELTA (desplazamiento DF) | 0.01 Å |
| Tamaños de supercelda DF | 2×2×2, 3×3×3, 4×4×4, 6×6×6 |
| k-grid de la celda primitiva | 12×12×12 |
| Malla q DFPT | 6×6×6 |
| Camino q | W → X → Γ → L → K |

El k-grid del SCF de la supercelda escala como `12 / N` (12×12×12 para
N=1, 2×2×2 para N=6), de modo que la densidad efectiva de muestreo en la
celda primitiva sea la misma en todos los tamaños.

## Resultados — modo óptico TO en Γ (×3 degenerado, cm⁻¹)

| N | DF | DFPT | Error relativo |
|---|-----|------|----------------|
| 2 | 473.5 | 467.01 | +1.4 % |
| 3 | 468.1 | 467.01 | +0.2 % |
| 4 | 467.3 | 467.01 | +0.06 % |
| 6 | **466.96** | **467.01** | **−0.01 %** |

DF converge suavemente a DFPT. En la mayor supercelda los dos métodos
coinciden por debajo del 0.01 %. La comparación de bandas completa está
en `PLOTS/DF_vs_DFPT_per_size.pdf`.

## Notas específicas

- Si es el caso "limpio" del proyecto: aislante, sin smearing, sin
  ambigüedades de funcional. **Usa LDA tanto en DF como en DFPT** (los
  dos pseudopotenciales son el mismo archivo, `pseudopotentzialak/Si.upf`).
- El cálculo DF en 6×6×6 (216 átomos primitivos × 2 = 432 atomos en la
  supercelda) es el más costoso de los tres materiales en términos
  computacionales (véase `denborak.csv`).
- No se aplica reducción por simetría: hay 12 SCFs desplazados por
  supercelda (6 direcciones × 2 átomos de la base).

## Reproducir

```bash
# Reconstruir bandas DF desde los outputs SCF ya incluidos (no requiere QE)
cd DF/6x6x6/out
python IFC_murriztue.py
python D.py
bash konparaketa.sh

# Agregar todos los tamaños
cd ../../../PLOTS
python datuak_batu.py
python plot_all.py
python erroreak.py

# Tiempos
cd ..
python denborak.py
python grafikoa_denb.py
```
