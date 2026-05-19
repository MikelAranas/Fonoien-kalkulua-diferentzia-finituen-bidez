#!/bin/bash
# =============================================================================
# convergence_ecutwfc.sh
# Test de convergencia de ecutwfc para Quantum ESPRESSO (pw.x)
#
# Criterio: comparación entre pasos CONSECUTIVOS (de menor a mayor ecutwfc)
#           |E(n) - E(n-1)| < THRESHOLD_MRY  durante DOS pasos seguidos
#
# Umbral: 0.1 mRy
# =============================================================================

# ---------------------------------------------------------------------------
# CONFIGURACIÓN
# ---------------------------------------------------------------------------
PW_EXEC="mpirun -np 12 pw.x"
PREFIX="graphene"
OUTDIR="./tmp"
PSEUDO_DIR="./"
PSEUDO="C.upf"

# Barrido de menor a mayor (la energía converge monótonamente al subir)
ECUT_LIST="20 30 40 50 60 80 100 120 140 160 180 200"

# Umbral en mRy — se para cuando DOS pasos consecutivos están por debajo
THRESHOLD_MRY=0.1

RESULTS_FILE="convergence_results.dat"
LOG_FILE="convergence_log.txt"

# ---------------------------------------------------------------------------
# Colores → siempre a stderr para no contaminar capturas de stdout
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# write_input ecut  — escribe pw_ecut<N>.in, sin output
# ---------------------------------------------------------------------------
write_input() {
    local ecut=$1
    cat > "pw_ecut${ecut}.in" << EOF
&CONTROL
  calculation = 'scf',
  prefix      = '${PREFIX}',
  outdir      = '${OUTDIR}',
  pseudo_dir  = '${PSEUDO_DIR}',
/
&SYSTEM
  ibrav     = 4,
  a         = 2.465325374,
  c         = 17.80,
  nat       = 2,
  ntyp      = 1,
  ecutwfc   = ${ecut}.0,
  occupations = 'smearing',
  smearing    = 'mv',
  degauss     = 0.02,
  assume_isolated = '2D'
/
&ELECTRONS
  mixing_beta = 0.7
  conv_thr    = 1.0d-9
/
ATOMIC_SPECIES
  C  12.0107 ${PSEUDO}
ATOMIC_POSITIONS (crystal)
  C    0.616642722    0.283357278   0.500000
  C    0.283357278    0.616642722   0.500000
K_POINTS automatic
  28 28 1   0 0 0
EOF
}

# ---------------------------------------------------------------------------
# run_pw ecut
#   stdout → solo la energía (lo que captura $())
#   stderr → mensajes de progreso (van a pantalla, no se capturan)
#   exit 0 → ok | exit 1 → error
# ---------------------------------------------------------------------------
run_pw() {
    local ecut=$1
    local infile="pw_ecut${ecut}.in"
    local outfile="pw_ecut${ecut}.out"

    echo -e "${CYAN}--- Calculando ecutwfc = ${ecut} Ry ---${NC}" >&2

    write_input "$ecut"

    echo -e "\n=== ecutwfc = ${ecut} Ry ===" >> "$LOG_FILE"
    > "$outfile"
    $PW_EXEC < "$infile" >> "$outfile" 2>&1
    local exit_code=$?
    cat "$outfile" >> "$LOG_FILE"

    if [ $exit_code -ne 0 ]; then
        echo -e "  ${RED}ERROR: pw.x falló (exit code ${exit_code}). Revisa: ${outfile}${NC}" >&2
        return 1
    fi

    local energy
    energy=$(grep "!    total energy" "$outfile" | tail -1 | awk '{print $5}')

    if [ -z "$energy" ]; then
        echo -e "  ${RED}AVISO: no se encontró 'total energy' en ${outfile}${NC}" >&2
        return 1
    fi

    echo "$energy"   # único output a stdout
    return 0
}

# ---------------------------------------------------------------------------
# delta_mry E1 E2  →  |E1-E2| en mRy  (stdout)
# ---------------------------------------------------------------------------
delta_mry() {
    awk -v a="$1" -v b="$2" \
        'BEGIN { d=a-b; if(d<0) d=-d; printf "%.5f", d*1000 }'
}

# ---------------------------------------------------------------------------
# is_below_threshold dmry  →  "SI" o "NO"
# ---------------------------------------------------------------------------
is_below() {
    awk -v d="$1" -v t="$THRESHOLD_MRY" \
        'BEGIN { print (d < t) ? "SI" : "NO" }'
}

# ---------------------------------------------------------------------------
# INICIO
# ---------------------------------------------------------------------------
echo -e "${CYAN}${BOLD}" >&2
echo "================================================================" >&2
echo "   Test de convergencia: ecutwfc — Quantum ESPRESSO"            >&2
echo "   Método  : comparación entre pasos CONSECUTIVOS"              >&2
echo "   Dirección: de menor a mayor ecutwfc"                         >&2
echo "   Umbral  : |E(n) - E(n-1)| < ${THRESHOLD_MRY} mRy"           >&2
echo "   Parada  : CINCO pasos consecutivos bajo el umbral"             >&2
echo "   conv_thr SCF: 1×10⁻⁹ Ry (un orden bajo el umbral)"          >&2
echo "================================================================" >&2
echo -e "${NC}" >&2

mkdir -p "$OUTDIR"
> "$LOG_FILE"

cat > "$RESULTS_FILE" << EOF
# Convergencia ecutwfc — pasos consecutivos   umbral: ${THRESHOLD_MRY} mRy
# conv_thr SCF: 1e-9 Ry
#
# ecutwfc(Ry)   E_total(Ry)            |ΔE|(mRy)       <umbral?   consec_ok
EOF

# Cabecera tabla pantalla
printf "\n${BOLD}%-14s %-24s %-14s %-12s %s${NC}\n" \
    "ecutwfc (Ry)" "E_total (Ry)" "|ΔE| (mRy)" "< umbral?" "Consec. ok" >&2
echo "------------------------------------------------------------------------" >&2

# ---------------------------------------------------------------------------
# BARRIDO
# ---------------------------------------------------------------------------
prev_energy=""
prev_below="NO"       # ¿el paso anterior estaba bajo el umbral?
consec=0              # contador de pasos consecutivos bajo el umbral
converged_at=""
streak_start=""

for ecut in $ECUT_LIST; do

    energy=$(run_pw "$ecut")
    if [ $? -ne 0 ] || [ -z "$energy" ]; then
        printf "%-14s %-24s %-14s %-12s %s\n" "$ecut" "ERROR" "---" "---" "---" >&2
        printf "%-14s %-24s %-14s %-12s %s\n" "$ecut" "ERROR" "---" "---" "---" \
            >> "$RESULTS_FILE"
        prev_energy=""
        prev_below="NO"
        consec=0
        continue
    fi

    if [ -z "$prev_energy" ]; then
        # Primer punto — solo mostramos energía, sin delta
        printf "%-14s %-24s %-14s %-12s %s\n" \
            "$ecut" "$energy" "---" "---" "(primero)" >&2
        printf "%-14s %-24s %-14s %-12s %s\n" \
            "$ecut" "$energy" "---" "---" "primero" >> "$RESULTS_FILE"
        prev_energy="$energy"
        continue
    fi

    dmry=$(delta_mry "$energy" "$prev_energy")
    below=$(is_below "$dmry")

    if [ "$below" = "SI" ]; then
        consec=$((consec + 1))
        if [ $consec -eq 1 ]; then
            streak_start=$ecut   # primer paso de la racha actual
        fi
    else
        consec=0
        streak_start=""
    fi

    # Etiqueta de pasos consecutivos
    if [ "$below" = "SI" ]; then
        consec_label="${GREEN}${consec}/5${NC}"
        consec_file="${consec}/5"
    else
        consec_label="${RED}0/5${NC}"
        consec_file="0/5"
    fi

    # Línea de tabla
    if [ "$below" = "SI" ]; then
        printf "%-14s %-24s %-14s " "$ecut" "$energy" "$dmry" >&2
        printf "${GREEN}%-12s${NC} " "✓ SI" >&2
        echo -e "$consec_label" >&2
    else
        printf "%-14s %-24s %-14s " "$ecut" "$energy" "$dmry" >&2
        printf "${RED}%-12s${NC} " "✗ NO" >&2
        echo -e "$consec_label" >&2
    fi

    printf "%-14s %-24s %-14s %-12s %s\n" \
        "$ecut" "$energy" "$dmry" "$below" "$consec_file" >> "$RESULTS_FILE"

    # Criterio de parada: 5 consecutivos bajo el umbral
    if [ $consec -ge 5 ] && [ -z "$converged_at" ]; then
        converged_at=$streak_start   # primero de los cinco pasos convergidos
        echo -e "\n  ${GREEN}${BOLD}✓ Convergencia alcanzada en ecutwfc = ${converged_at} Ry${NC}" >&2
        echo -e "  (cinco pasos consecutivos con |ΔE| < ${THRESHOLD_MRY} mRy)" >&2
        echo -e "  Continuando para confirmar estabilidad...\n" >&2
    fi

    prev_energy="$energy"
    prev_ecut="$ecut"
    prev_below="$below"

done

# ---------------------------------------------------------------------------
# RESUMEN FINAL
# ---------------------------------------------------------------------------
echo "" >&2
echo "================================================================" >&2
echo -e "${BOLD}RESUMEN FINAL${NC}" >&2
echo "----------------------------------------------------------------" >&2

if [ -n "$converged_at" ]; then
    echo -e "  ${GREEN}${BOLD}Convergencia alcanzada: ecutwfc = ${converged_at} Ry${NC}" >&2
    echo -e "  Criterio: |ΔE| < ${THRESHOLD_MRY} mRy en cinco pasos consecutivos" >&2
    echo "" >&2
    echo -e "  → Usa ${GREEN}${BOLD}ecutwfc = ${converged_at} Ry${NC} en tus cálculos." >&2
else
    echo -e "  ${RED}No se alcanzó convergencia con los valores probados.${NC}" >&2
    echo -e "  Añade valores más altos a ECUT_LIST." >&2
fi

echo "" >&2
echo "  Tabla de resultados : $RESULTS_FILE" >&2
echo "  Log completo pw.x   : $LOG_FILE"     >&2
echo "================================================================" >&2
