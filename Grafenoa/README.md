# Grafeno

Semimetal 2D, red en panal de abeja (hexagonal con 2 átomos en la base).
**Caso más exigente del proyecto** — tanto por la naturaleza 2D como por
la sensibilidad de los modos ópticos al funcional XC.

## Parámetros

| Magnitud | Valor |
|----------|-------|
| Constante de red en el plano `a` | 2.4653 Å |
| Vacío fuera del plano `c` | 17.80 Å |
| Base | C en (0, 0, 0), C en (1/3, 2/3, 0) |
| Masa atómica | 12.011 uma |
| Pseudopotencial | `C.upf` — ATOMPAW, **LDA** |
| Cutoff (ecutwfc) | 100 Ry |
| Ocupaciones | `smearing = 'mv'` (Marzari–Vanderbilt), `degauss = 0.02` Ry |
| `assume_isolated` | `2D` (truncación de Coulomb en z) |
| `conv_thr` (SCF) | 10⁻¹⁶ Ry |
| DELTA (desplazamiento DF) | 0.01 Å |
| Tamaños de supercelda DF | 2×2×1 (8 at), 4×4×1 (32 at), 7×7×1 (98 at) |
| k-grid de la celda primitiva | 28×28×1 |
| Malla q DFPT | 7×7×1 |
| Camino q | Γ → M → K → Γ |

**Detalles importantes**:

- `assume_isolated = '2D'` es obligatorio en QE para que la interacción
  electrostática entre imágenes periódicas en la dirección z no
  contamine los cálculos.
- Nunca se usan k-grids ni mallas q múltiplos de 3: caerían en el punto
  K (Dirac), donde la estructura electrónica es singular y produce
  inestabilidades fonónicas.
- En `D.py`, las posiciones atómicas de la base se construyen como
  **fracciones exactas de Python** (`1/3`, `2/3`), no decimales
  truncados. Una truncación de 10⁻⁶ basta para que la multiplicidad de
  un átomo en la frontera de Wigner–Seitz baje de 3 a 2, lo que rompe la
  degeneración del cono de Dirac en K (separación de ~26 cm⁻¹).

## Resultados — frecuencias en Γ (cm⁻¹)

DFPT y DF usan **ambos el pseudopotencial LDA** (un único archivo,
`pseudopotentzialak/C.upf`).

| Supercelda | ZO (DF) | ZO (DFPT) | LO/TO (DF) | LO/TO (DFPT) |
|------------|---------|-----------|------------|--------------|
| 2×2×1 | 900.94 | 906.67 | 1545.82 / 1545.84 | 1545.60 / 1545.60 |
| 4×4×1 | 901.00 | 902.52 | 1545.95 / 1545.97 | 1545.04 / 1545.04 |
| 7×7×1 | **901.90** | **906.67** | **1546.13 / 1546.16** | **1545.60 / 1545.60** |

Las bandas completas están en `PLOTS/DF_vs_DFPT_grafeno.pdf` y el
análisis de error en `PLOTS/errore_erlatiboa_grafeno.pdf`.

## Limitaciones conocidas

- **Rama ZA en DF**: el modo acústico fuera del plano cerca de Γ sale
  *lineal* con DF, no cuadrático como pide la física del grafeno
  (ω ∝ q²). Es un artefacto intrínseco de la interpolación 3D-periódica
  del IFC en espacio real. Para comparar la rama ZA hay que usar DFPT
  con malla q suficientemente fina (7×7×1 funciona, 4×4×1 no).
- **Rama ZA negativa entre puntos de malla en DFPT 4×4×1**: matdyn
  interpola entre los 16 puntos q de la malla y para el modo flexural
  (q²) la interpolación se vuelve inestable, generando frecuencias
  imaginarias en `q < 0.3 |G|` a lo largo de Γ → M. **No es física**:
  los valores DFPT *en los puntos de la malla* son todos positivos. Se
  resuelve usando una malla q más densa.

## Lecciones metodológicas (resumen)

Estas dos lecciones costaron buena parte del debugging del proyecto:

1. **DF y DFPT deben usar el mismo `.upf`.** El nombre del archivo no
   garantiza el contenido. Una semana se perdió atribuyendo un
   desacuerdo de 28 cm⁻¹ a una "limitación intrínseca del método DF en
   2D" cuando en realidad DF usaba LDA y DFPT usaba PBE. El `md5sum` y
   la línea `Exchange-correlation` del output de `pw.x` son las
   verificaciones definitivas.
2. **Las fracciones de la base deben ser exactas en `D.py`.** El átomo
   B del grafeno está en (1/3, 2/3); escribirlo como `0.333333` rompe
   la multiplicidad en la frontera de Wigner–Seitz y se pierde la
   degeneración del cono de Dirac en K.

## Reproducir

```bash
cd DF/7x7x1/out
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
