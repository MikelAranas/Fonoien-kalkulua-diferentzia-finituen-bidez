#!/bin/bash
# =============================================================================
# kgrid.sh — Test de convergencia de k-grid para Si (pw.x)
# Adaptado de Grafenoa/DFPT/conv/kgrid.sh
#
# Criterion: comparison between CONSECUTIVE grids (low to high)
#            |E(n) - E(n-1)| < THRESHOLD_MRY for FIVE steps in a row
#
# ⚠ ECUTWFC fijo al valor convergido por ecut.sh — EDITAR ANTES DE LANZAR.
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

ECUTWFC=50.0          

# Grids to test, from low to high density (malla N N N sin shift)
KGRID_LIST="2 4 6 8 10 12 14 16 18 20"

THRESHOLD_MRY=0.1

RESULTS_FILE="convergence_kpoints.dat"
LOG_FILE="convergence_kpoints_log.txt"

RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

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
  ibrav     = 2,
  a         = ${ALAT},
  nat       = 2,
  ntyp      = 1,
  ecutwfc   = ${ECUTWFC},
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
  ${n} ${n} ${n}   0 0 0
EOF
}

run_pw() {
    local n=$1
    local infile="pw_k${n}.in"
    local outfile="pw_k${n}.out"

    echo -e "${CYAN}--- Calculando k-grid = ${n}x${n}x${n} ---${NC}" >&2

    write_input "$n"

    echo -e "\n=== k-grid = ${n}x${n}x${n} ===" >> "$LOG_FILE"
    > "$outfile"
    $PW_EXEC < "$infile" >> "$outfile" 2>&1
    local exit_code=$?
    cat "$outfile" >> "$LOG_FILE"

    if [ $exit_code -ne 0 ]; then
        echo -e "  ${RED}ERROR: pw.x falló (exit ${exit_code}). Revisa: ${outfile}${NC}" >&2
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
echo -e "${CYAN}${BOLD}" >&2
echo "================================================================" >&2
echo "   Test de convergencia: k-grid — Si (NC-PBE)"                   >&2
echo "   ecutwfc fijo: ${ECUTWFC} Ry"                                  >&2
echo "   Threshold: |E(n) - E(n-1)| < ${THRESHOLD_MRY} mRy"             >&2
echo "   Stop    : FIVE consecutive steps below threshold"            >&2
echo "================================================================" >&2
echo -e "${NC}" >&2

mkdir -p "$OUTDIR"
> "$LOG_FILE"

cat > "$RESULTS_FILE" << EOF
# Convergencia k-grid — ecutwfc: ${ECUTWFC} Ry   umbral: ${THRESHOLD_MRY} mRy
# Si FCC diamante, a=${ALAT} A, conv_thr SCF: 1e-9 Ry
#
# k-grid       E_total(Ry)            |ΔE|(mRy)       <umbral?   consec
EOF

printf "\n${BOLD}%-12s %-24s %-14s %-12s %s${NC}\n" \
    "k-grid" "E_total (Ry)" "|ΔE| (mRy)" "< umbral?" "Consec. ok" >&2
echo "------------------------------------------------------------------------" >&2

prev_energy=""
consec=0
converged_at=""
streak_start=""

for n in $KGRID_LIST; do

    label="${n}x${n}x${n}"
    energy=$(run_pw "$n")
    if [ $? -ne 0 ] || [ -z "$energy" ]; then
        printf "%-12s %-24s %-14s %-12s %s\n" "$label" "ERROR" "---" "---" "---" >&2
        printf "%-12s %-24s %-14s %-12s %s\n" "$label" "ERROR" "---" "---" "---" \
            >> "$RESULTS_FILE"
        prev_energy=""
        consec=0
        continue
    fi

    if [ -z "$prev_energy" ]; then
        printf "%-12s %-24s %-14s %-12s %s\n" \
            "$label" "$energy" "---" "---" "(primero)" >&2
        printf "%-12s %-24s %-14s %-12s %s\n" \
            "$label" "$energy" "---" "---" "primero" >> "$RESULTS_FILE"
        prev_energy="$energy"
        continue
    fi

    dmry=$(delta_mry "$energy" "$prev_energy")
    below=$(is_below "$dmry")

    if [ "$below" = "SI" ]; then
        consec=$((consec + 1))
        if [ $consec -eq 1 ]; then
            streak_start=$label
        fi
    else
        consec=0
        streak_start=""
    fi

    if [ "$below" = "SI" ]; then
        printf "%-12s %-24s %-14s ${GREEN}%-12s${NC} ${GREEN}%s/5${NC}\n" \
            "$label" "$energy" "$dmry" "✓ SI" "$consec" >&2
        consec_file="${consec}/5"
    else
        printf "%-12s %-24s %-14s ${RED}%-12s${NC} ${RED}0/5${NC}\n" \
            "$label" "$energy" "$dmry" "✗ NO" >&2
        consec_file="0/5"
    fi

    printf "%-12s %-24s %-14s %-12s %s\n" \
        "$label" "$energy" "$dmry" "$below" "$consec_file" >> "$RESULTS_FILE"

    if [ $consec -ge 5 ] && [ -z "$converged_at" ]; then
        converged_at=$streak_start
        echo -e "\n  ${GREEN}${BOLD}✓ Convergencia alcanzada en k-grid = ${converged_at}${NC}" >&2
        echo -e "  Continuing to confirm stability...\n" >&2
    fi

    prev_energy="$energy"

done

echo "" >&2
echo "================================================================" >&2
echo -e "${BOLD}FINAL SUMMARY${NC}" >&2
echo "----------------------------------------------------------------" >&2
if [ -n "$converged_at" ]; then
    echo -e "  ${GREEN}${BOLD}Convergence reached: k-grid = ${converged_at}${NC}" >&2
    echo -e "  → Edita KGRID en relax.sh y úsalo en si.scf.in" >&2
else
    echo -e "  ${RED}Convergence not reached with the tested values.${NC}" >&2
fi
echo "" >&2
echo "  Results table     : $RESULTS_FILE" >&2
echo "  Full pw.x log     : $LOG_FILE"     >&2
echo "================================================================" >&2
