#!/bin/bash
# =============================================================================
# convergence_ecutwfc.sh
# ecutwfc convergence test for Quantum ESPRESSO (pw.x)
#
# Criterion: comparison between CONSECUTIVE steps (low to high ecutwfc)
#           |E(n) - E(n-1)| < THRESHOLD_MRY  for TWO steps in a row
#
# Threshold: 0.1 mRy
# =============================================================================

# ---------------------------------------------------------------------------
# CONFIGURATION
# ---------------------------------------------------------------------------
PW_EXEC="mpirun -np 48 pw.x"
PREFIX="al"
OUTDIR="./tmp"
PSEUDO_DIR="../../../pseudopotentzialak/"
PSEUDO="Al.upf"

# Sweep low to high (energy converges monotonically as it increases)
ECUT_LIST="20 30 40 50 60 80 100 120 140 160 180 200"

# Threshold in mRy - stops when TWO consecutive steps are below it
THRESHOLD_MRY=0.1

RESULTS_FILE="convergence_results.dat"
LOG_FILE="convergence_log.txt"

# ---------------------------------------------------------------------------
# Colors -> always to stderr so stdout captures stay clean
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
YELLOW='\033[1;33m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# write_input ecut  — writes pw_ecut<N>.in, no output
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
  ibrav     = 0,
  nat       = 1,
  ntyp      = 1,
  ecutwfc   = ${ecut}.0,
  occupations = 'smearing',
  smearing    = 'mp',
  degauss     = 0.02,
/
&ELECTRONS
  mixing_beta = 0.7
  conv_thr    = 1.0d-16
/
ATOMIC_SPECIES
  Al  26.9815 ${PSEUDO}
ATOMIC_POSITIONS angstrom
  Al  0.00 0.00 0.00
K_POINTS automatic
  12 12 12 0 0 0
CELL_PARAMETERS angstrom
  -2.0189  0.0000  2.0189
   0.0000  2.0189  2.0189
  -2.0189  2.0189  0.0000
EOF
}

# ---------------------------------------------------------------------------
# run_pw ecut
#   stdout -> only the energy (what $() captures)
#   stderr -> progress messages (go to screen, not captured)
#   exit 0 -> ok | exit 1 -> error
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

    echo "$energy"   # only output to stdout
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
echo "   Method  : comparison between CONSECUTIVE steps"              >&2
echo "   Direction: low to high ecutwfc"                         >&2
echo "   Threshold: |E(n) - E(n-1)| < ${THRESHOLD_MRY} mRy"           >&2
echo "   Stop    : FIVE consecutive steps below threshold"             >&2
echo "   SCF conv_thr: 1e-9 Ry (one order below the threshold)"          >&2
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

# On-screen table header
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
        # First point - only show energy, no delta
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
            streak_start=$ecut   # first step of the current streak
        fi
    else
        consec=0
        streak_start=""
    fi

    # Consecutive-steps label
    if [ "$below" = "SI" ]; then
        consec_label="${GREEN}${consec}/5${NC}"
        consec_file="${consec}/5"
    else
        consec_label="${RED}0/5${NC}"
        consec_file="0/5"
    fi

    # Table row
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

    # Stop criterion: 5 consecutive below threshold
    if [ $consec -ge 5 ] && [ -z "$converged_at" ]; then
        converged_at=$streak_start   # first of the five converged steps
        echo -e "\n  ${GREEN}${BOLD}✓ Convergence reached at ecutwfc = ${converged_at} Ry${NC}" >&2
        echo -e "  (five consecutive steps with |ΔE| < ${THRESHOLD_MRY} mRy)" >&2
        echo -e "  Continuing to confirm stability...\n" >&2
    fi

    prev_energy="$energy"
    prev_ecut="$ecut"
    prev_below="$below"

done

# ---------------------------------------------------------------------------
# FINAL SUMMARY
# ---------------------------------------------------------------------------
echo "" >&2
echo "================================================================" >&2
echo -e "${BOLD}FINAL SUMMARY${NC}" >&2
echo "----------------------------------------------------------------" >&2

if [ -n "$converged_at" ]; then
    echo -e "  ${GREEN}${BOLD}Convergence reached: ecutwfc = ${converged_at} Ry${NC}" >&2
    echo -e "  Criterion: |dE| < ${THRESHOLD_MRY} mRy over five consecutive steps" >&2
    echo "" >&2
    echo -e "  → Usa ${GREEN}${BOLD}ecutwfc = ${converged_at} Ry${NC} in your calculations." >&2
else
    echo -e "  ${RED}Convergence not reached with the tested values.${NC}" >&2
    echo -e "  Add higher values to ECUT_LIST." >&2
fi

echo "" >&2
echo "  Results table     : $RESULTS_FILE" >&2
echo "  Full pw.x log     : $LOG_FILE"     >&2
echo "================================================================" >&2
