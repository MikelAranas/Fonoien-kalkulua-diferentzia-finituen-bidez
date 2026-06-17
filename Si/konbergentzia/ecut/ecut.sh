#!/bin/bash
# =============================================================================
# ecut.sh — Test de convergencia de ecutwfc para Si (pw.x)
# Adaptado de Grafenoa/DFPT/conv/ecut.sh
#
# Criterion: comparison between CONSECUTIVE steps (low to high ecutwfc)
#           |E(n) - E(n-1)| < THRESHOLD_MRY for FIVE steps in a row
# Threshold: 0.1 mRy
#
# Geometría fija: FCC diamante, a = 5.431 Å (experimental; el a definitivo
# sale después de relax.sh). k-grid fijo 12x12x12 (densidad de producción).
# =============================================================================

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
PW_EXEC="mpirun -np 48 pw.x"
PREFIX="si"
OUTDIR="./tmp"
PSEUDO_DIR="../../../pseudopotentzialak/"
PSEUDO="Si_NC_PBE.upf"
ALAT=5.431

# Sweep low to high (energy converges monotonically as it increases)
ECUT_LIST="15 20 25 30 35 40 45 50 60 70 80 90 100 110 120 130 140 150 160"

THRESHOLD_MRY=0.1

RESULTS_FILE="convergence_results.dat"
LOG_FILE="convergence_log.txt"

# ---------------------------------------------------------------------------
# Colors -> always to stderr so stdout captures stay clean
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# write_input ecut — writes pw_ecut<N>.in, no output
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
  ibrav     = 2,
  a         = ${ALAT},
  nat       = 2,
  ntyp      = 1,
  ecutwfc   = ${ecut}.0,
  occupations = 'fixed'
/
&ELECTRONS
  conv_thr    = 1.0d-8
/
ATOMIC_SPECIES
  Si  28.086 ${PSEUDO}
ATOMIC_POSITIONS (crystal)
  Si   0.000000   0.000000   0.000000
  Si   0.250000   0.250000   0.250000
K_POINTS automatic
  12 12 12   0 0 0
EOF
}

# ---------------------------------------------------------------------------
# run_pw ecut
#   stdout → solo la energía | stderr → progreso
# ---------------------------------------------------------------------------
run_pw() {
    local ecut=$1
    local infile="pw_ecut${ecut}.in"
    local outfile="pw_ecut${ecut}.out"

    echo -e "${CYAN}--- Computing ecutwfc = ${ecut} Ry ---${NC}" >&2

    write_input "$ecut"

    echo -e "\n=== ecutwfc = ${ecut} Ry ===" >> "$LOG_FILE"
    > "$outfile"
    $PW_EXEC < "$infile" >> "$outfile" 2>&1
    local exit_code=$?
    cat "$outfile" >> "$LOG_FILE"

    if [ $exit_code -ne 0 ]; then
        echo -e "  ${RED}ERROR: pw.x failed (exit code ${exit_code}). Revisa: ${outfile}${NC}" >&2
        return 1
    fi

    local energy
    energy=$(grep "!    total energy" "$outfile" | tail -1 | awk '{print $5}')

    if [ -z "$energy" ]; then
        echo -e "  ${RED}WARNING: 'total energy' not found in ${outfile}${NC}" >&2
        return 1
    fi

    echo "$energy"
    return 0
}

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
echo "   Test de convergencia: ecutwfc — Si (NC-PBE)"                  >&2
echo "   Method  : comparison between CONSECUTIVE steps"               >&2
echo "   Threshold: |E(n) - E(n-1)| < ${THRESHOLD_MRY} mRy"             >&2
echo "   Stop    : FIVE consecutive steps below threshold"            >&2
echo "   conv_thr SCF: 1×10⁻⁹ Ry"                                      >&2
echo "================================================================" >&2
echo -e "${NC}" >&2

mkdir -p "$OUTDIR"
> "$LOG_FILE"

cat > "$RESULTS_FILE" << EOF
# Convergencia ecutwfc — pasos consecutivos   umbral: ${THRESHOLD_MRY} mRy
# Si FCC diamante, a=${ALAT} A, k=12x12x12, conv_thr SCF: 1e-9 Ry
#
# ecutwfc(Ry)   E_total(Ry)            |ΔE|(mRy)       <umbral?   consec_ok
EOF

printf "\n${BOLD}%-14s %-24s %-14s %-12s %s${NC}\n" \
    "ecutwfc (Ry)" "E_total (Ry)" "|ΔE| (mRy)" "< umbral?" "Consec. ok" >&2
echo "------------------------------------------------------------------------" >&2

prev_energy=""
consec=0
converged_at=""
streak_start=""

for ecut in $ECUT_LIST; do

    energy=$(run_pw "$ecut")
    if [ $? -ne 0 ] || [ -z "$energy" ]; then
        printf "%-14s %-24s %-14s %-12s %s\n" "$ecut" "ERROR" "---" "---" "---" >&2
        printf "%-14s %-24s %-14s %-12s %s\n" "$ecut" "ERROR" "---" "---" "---" \
            >> "$RESULTS_FILE"
        prev_energy=""
        consec=0
        continue
    fi

    if [ -z "$prev_energy" ]; then
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
            streak_start=$ecut
        fi
    else
        consec=0
        streak_start=""
    fi

    if [ "$below" = "SI" ]; then
        printf "%-14s %-24s %-14s ${GREEN}%-12s${NC} ${GREEN}%s/5${NC}\n" \
            "$ecut" "$energy" "$dmry" "✓ SI" "$consec" >&2
        consec_file="${consec}/5"
    else
        printf "%-14s %-24s %-14s ${RED}%-12s${NC} ${RED}0/5${NC}\n" \
            "$ecut" "$energy" "$dmry" "✗ NO" >&2
        consec_file="0/5"
    fi

    printf "%-14s %-24s %-14s %-12s %s\n" \
        "$ecut" "$energy" "$dmry" "$below" "$consec_file" >> "$RESULTS_FILE"

    if [ $consec -ge 5 ] && [ -z "$converged_at" ]; then
        converged_at=$streak_start
        echo -e "\n  ${GREEN}${BOLD}✓ Convergence reached at ecutwfc = ${converged_at} Ry${NC}" >&2
        echo -e "  Continuing to confirm stability...\n" >&2
    fi

    prev_energy="$energy"

done

echo "" >&2
echo "================================================================" >&2
echo -e "${BOLD}FINAL SUMMARY${NC}" >&2
echo "----------------------------------------------------------------" >&2
if [ -n "$converged_at" ]; then
    echo -e "  ${GREEN}${BOLD}Convergence reached: ecutwfc = ${converged_at} Ry${NC}" >&2
    echo -e "  → Edita ECUTWFC=${converged_at}.0 en kgrid.sh y relax.sh" >&2
else
    echo -e "  ${RED}Convergence not reached with the tested values.${NC}" >&2
    echo -e "  Add higher values to ECUT_LIST." >&2
fi
echo "" >&2
echo "  Results table     : $RESULTS_FILE" >&2
echo "  Full pw.x log     : $LOG_FILE"     >&2
echo "================================================================" >&2
