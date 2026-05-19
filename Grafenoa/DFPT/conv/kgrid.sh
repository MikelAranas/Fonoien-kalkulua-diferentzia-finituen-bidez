#!/bin/bash
# =============================================================================
# convergence_kpoints.sh
# Test de convergencia de k-grid para Quantum ESPRESSO (pw.x) — grafeno 2D
#
# Criterio : comparación entre grids CONSECUTIVOS (de menor a mayor)
#            |E(n) - E(n-1)| < THRESHOLD_MRY  durante CINCO pasos seguidos
#
# ecutwfc fijo al valor convergido
# =============================================================================

# ---------------------------------------------------------------------------
# CONFIGURACIÓN
# ---------------------------------------------------------------------------
PW_EXEC="mpirun -np 12 pw.x"
PREFIX="graphene"
OUTDIR="./tmp"
PSEUDO_DIR="./"
PSEUDO="C.upf"

ECUTWFC=100.0          # valor convergido del test anterior

# Grids a probar, de menor a mayor densidad
# Formato: "NxN" — el script construye "N N 1  0 0 0" automáticamente
KGRID_LIST="4 8 12 16 20 24 28 32 36 40 44 48"

THRESHOLD_MRY=0.1
NSTEPS=5               # pasos consecutivos requeridos

RESULTS_FILE="convergence_kpoints.dat"
LOG_FILE="convergence_kpoints_log.txt"

# ---------------------------------------------------------------------------
# Colores → stderr para no contaminar capturas de stdout
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# write_input N  — escribe pw_kN.in
# ---------------------------------------------------------------------------
write_input() {
    local n=$1
    cat > "pw_k${n}.in" << EOF
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
  ecutwfc   = ${ECUTWFC},
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
  ${n} ${n} 1   0 0 0
EOF
}

# ---------------------------------------------------------------------------
# run_pw N
#   stdout → solo la energía
#   stderr → mensajes de progreso
# ---------------------------------------------------------------------------
run_pw() {
    local n=$1
    local infile="pw_k${n}.in"
    local outfile="pw_k${n}.out"

    echo -e "${CYAN}--- Calculando k-grid = ${n}x${n}x1 ---${NC}" >&2

    write_input "$n"

    echo -e "\n=== k-grid = ${n}x${n}x1 ===" >> "$LOG_FILE"
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

    echo "$energy"
    return 0
}

# ---------------------------------------------------------------------------
# delta_mry E1 E2  →  |E1-E2| en mRy
# ---------------------------------------------------------------------------
delta_mry() {
    awk -v a="$1" -v b="$2" \
        'BEGIN { d=a-b; if(d<0) d=-d; printf "%.5f", d*1000 }'
}

is_below() {
    awk -v d="$1" -v t="$THRESHOLD_MRY" \
        'BEGIN { print (d < t) ? "SI" : "NO" }'
}

# ---------------------------------------------------------------------------
# INICIO
# ---------------------------------------------------------------------------
echo -e "${CYAN}${BOLD}" >&2
echo "================================================================" >&2
echo "   Test de convergencia: k-grid — Quantum ESPRESSO"             >&2
echo "   ecutwfc fijo : ${ECUTWFC} Ry"                                >&2
echo "   Método       : comparación entre grids CONSECUTIVOS"         >&2
echo "   Dirección    : de menor a mayor densidad"                     >&2
echo "   Umbral       : |E(n) - E(n-1)| < ${THRESHOLD_MRY} mRy"      >&2
echo "   Parada       : ${NSTEPS} pasos consecutivos bajo el umbral"  >&2
echo "   conv_thr SCF : 1×10⁻⁹ Ry"                                   >&2
echo "================================================================" >&2
echo -e "${NC}" >&2

mkdir -p "$OUTDIR"
> "$LOG_FILE"

cat > "$RESULTS_FILE" << EOF
# Convergencia k-grid — ecutwfc: ${ECUTWFC} Ry   umbral: ${THRESHOLD_MRY} mRy
# conv_thr SCF: 1e-9 Ry
#
# k-grid       E_total(Ry)            |ΔE|(mRy)       <umbral?   consec
EOF

# Cabecera tabla pantalla
printf "\n${BOLD}%-14s %-24s %-14s %-12s %s${NC}\n" \
    "k-grid" "E_total (Ry)" "|ΔE| (mRy)" "< umbral?" "Consec. ok" >&2
echo "------------------------------------------------------------------------" >&2

# ---------------------------------------------------------------------------
# BARRIDO
# ---------------------------------------------------------------------------
prev_energy=""
consec=0
streak_start=""
converged_at=""

for n in $KGRID_LIST; do

    energy=$(run_pw "$n")
    if [ $? -ne 0 ] || [ -z "$energy" ]; then
        printf "%-14s %-24s %-14s %-12s %s\n" "${n}x${n}x1" "ERROR" "---" "---" "---" >&2
        printf "%-14s %-24s %-14s %-12s %s\n" "${n}x${n}x1" "ERROR" "---" "---" "---" \
            >> "$RESULTS_FILE"
        prev_energy=""
        consec=0
        streak_start=""
        continue
    fi

    if [ -z "$prev_energy" ]; then
        printf "%-14s %-24s %-14s %-12s %s\n" \
            "${n}x${n}x1" "$energy" "---" "---" "(primero)" >&2
        printf "%-14s %-24s %-14s %-12s %s\n" \
            "${n}x${n}x1" "$energy" "---" "---" "primero" >> "$RESULTS_FILE"
        prev_energy="$energy"
        continue
    fi

    dmry=$(delta_mry "$energy" "$prev_energy")
    below=$(is_below "$dmry")

    if [ "$below" = "SI" ]; then
        consec=$((consec + 1))
        if [ $consec -eq 1 ]; then
            streak_start="${n}x${n}x1"
        fi
    else
        consec=0
        streak_start=""
    fi

    if [ "$below" = "SI" ]; then
        consec_label="${GREEN}${consec}/${NSTEPS}${NC}"
        consec_file="${consec}/${NSTEPS}"
        printf "%-14s %-24s %-14s " "${n}x${n}x1" "$energy" "$dmry" >&2
        printf "${GREEN}%-12s${NC} " "✓ SI" >&2
        echo -e "$consec_label" >&2
    else
        consec_label="${RED}0/${NSTEPS}${NC}"
        consec_file="0/${NSTEPS}"
        printf "%-14s %-24s %-14s " "${n}x${n}x1" "$energy" "$dmry" >&2
        printf "${RED}%-12s${NC} " "✗ NO" >&2
        echo -e "$consec_label" >&2
    fi

    printf "%-14s %-24s %-14s %-12s %s\n" \
        "${n}x${n}x1" "$energy" "$dmry" "$below" "$consec_file" >> "$RESULTS_FILE"

    if [ $consec -ge $NSTEPS ] && [ -z "$converged_at" ]; then
        converged_at=$streak_start
        echo -e "\n  ${GREEN}${BOLD}✓ Convergencia alcanzada en k-grid = ${converged_at}${NC}" >&2
        echo -e "  (${NSTEPS} pasos consecutivos con |ΔE| < ${THRESHOLD_MRY} mRy)" >&2
        echo -e "  Continuando para confirmar estabilidad...\n" >&2
    fi

    prev_energy="$energy"

done

# ---------------------------------------------------------------------------
# RESUMEN FINAL
# ---------------------------------------------------------------------------
echo "" >&2
echo "================================================================" >&2
echo -e "${BOLD}RESUMEN FINAL${NC}" >&2
echo "----------------------------------------------------------------" >&2

if [ -n "$converged_at" ]; then
    echo -e "  ${GREEN}${BOLD}Convergencia alcanzada: k-grid = ${converged_at}${NC}" >&2
    echo -e "  Criterio: |ΔE| < ${THRESHOLD_MRY} mRy en ${NSTEPS} pasos consecutivos" >&2
    echo "" >&2
    echo -e "  → Usa ${GREEN}${BOLD}K_POINTS automatic  ${converged_at%x*} 1   0 0 0${NC}" >&2
else
    echo -e "  ${RED}No se alcanzó convergencia. Añade grids más densos a KGRID_LIST.${NC}" >&2
fi

echo "" >&2
echo "  Tabla de resultados : $RESULTS_FILE" >&2
echo "  Log completo pw.x   : $LOG_FILE"     >&2
echo "================================================================" >&2
