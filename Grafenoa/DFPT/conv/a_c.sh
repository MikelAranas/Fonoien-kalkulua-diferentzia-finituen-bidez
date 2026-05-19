#!/bin/bash
# =============================================================================
# lattice_sweep.sh
# Barrido de parámetro de red a y c para grafeno — Quantum ESPRESSO
# Fija uno mientras barre el otro, encuentra el mínimo de energía
# =============================================================================

# ---------------------------------------------------------------------------
# CONFIGURACIÓN
# ---------------------------------------------------------------------------
PW_EXEC="mpirun -np 12 pw.x"
PREFIX="graphene"
OUTDIR="./tmp"
PSEUDO_DIR="./"
PSEUDO="C.upf"

ECUTWFC=100.0
KGRID="28 28 1"

# Barrido de a (Å) — centro en valor experimental
A_CENTER=2.465325374
A_STEP=0.02
A_NSTEPS=6        # ±6 pasos → 2.345 a 2.585 Å

# Barrido de c (Å) — centro en valor actual
C_CENTER=17.80
C_STEP=0.5
C_NSTEPS=6        # ±6 pasos → 14.8 a 20.8 Å

# a fijo para el barrido de c (se actualiza tras encontrar mínimo en a)
A_FIXED=$A_CENTER
# c fijo para el barrido de a
C_FIXED=$C_CENTER

RESULTS_A="sweep_a.dat"
RESULTS_C="sweep_c.dat"
LOG_FILE="sweep_log.txt"

# ---------------------------------------------------------------------------
# Colores → stderr
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
BOLD='\033[1m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# write_input a c outfile
# ---------------------------------------------------------------------------
write_input() {
    local a=$1
    local c=$2
    local infile=$3

    # Vectores de red explícitos a partir de a y c
    local a1_x=$(awk -v a="$a" 'BEGIN{printf "%.10f", a}')
    local a2_x=$(awk -v a="$a" 'BEGIN{printf "%.10f", -a/2}')
    local a2_y=$(awk -v a="$a" 'BEGIN{printf "%.10f", a*sqrt(3)/2}')

    cat > "$infile" << EOF
&CONTROL
  calculation = 'scf',
  prefix      = '${PREFIX}',
  outdir      = '${OUTDIR}',
  pseudo_dir  = '${PSEUDO_DIR}',
/
&SYSTEM
  ibrav     = 0,
  nat       = 2,
  ntyp      = 1,
  ecutwfc   = ${ECUTWFC},
  occupations = 'smearing',
  smearing    = 'mv',
  degauss     = 0.02,
  assume_isolated = '2D'
/
&ELECTRONS
  mixing_beta = 0.7,
  conv_thr    = 1.0d-9,
/
ATOMIC_SPECIES
  C  12.0107  ${PSEUDO}

CELL_PARAMETERS (angstrom)
  ${a1_x}   0.0000000000   0.0000000000
  ${a2_x}   ${a2_y}   0.0000000000
  0.0000000000   0.0000000000   ${c}

ATOMIC_POSITIONS (crystal)
  C    0.616642722    0.283357278   0.500000
  C    0.283357278    0.616642722   0.500000

K_POINTS automatic
  ${KGRID}   0 0 0
EOF
}

# ---------------------------------------------------------------------------
# run_pw infile outfile → energía a stdout, mensajes a stderr
# ---------------------------------------------------------------------------
run_pw() {
    local infile=$1
    local outfile=$2

    $PW_EXEC < "$infile" > "$outfile" 2>&1
    local exit_code=$?
    cat "$outfile" >> "$LOG_FILE"

    if [ $exit_code -ne 0 ]; then
        echo -e "  ${RED}ERROR: pw.x falló. Revisa: ${outfile}${NC}" >&2
        return 1
    fi

    local energy
    energy=$(grep "!    total energy" "$outfile" | tail -1 | awk '{print $5}')
    if [ -z "$energy" ]; then
        echo -e "  ${RED}AVISO: no se encontró total energy en ${outfile}${NC}" >&2
        return 1
    fi

    echo "$energy"
    return 0
}

# ---------------------------------------------------------------------------
# find_min datfile → imprime el valor de param con menor energía
# ---------------------------------------------------------------------------
find_min() {
    # Formato del dat: param  energy
    awk 'NF==2 && $1!~/^#/ {
        if (NR==1 || $2 < min_e) { min_e=$2; min_p=$1 }
    } END { print min_p, min_e }' "$1"
}

# ---------------------------------------------------------------------------
# INICIO
# ---------------------------------------------------------------------------
echo -e "${CYAN}${BOLD}" >&2
echo "================================================================" >&2
echo "   Barrido de parámetros de red — Quantum ESPRESSO"             >&2
echo "   ecutwfc : ${ECUTWFC} Ry   k-grid : ${KGRID}"                >&2
echo "================================================================" >&2
echo -e "${NC}" >&2

mkdir -p "$OUTDIR"
> "$LOG_FILE"

# ===========================================================================
# PASO 1: barrido de a (c fijo)
# ===========================================================================
echo -e "\n${BOLD}[1/2] Barrido de a  (c fijo = ${C_FIXED} Å)${NC}" >&2
echo "# Barrido de a — c fijo = ${C_FIXED} A   ecutwfc=${ECUTWFC} Ry" > "$RESULTS_A"
echo "# a(Ang)      E_total(Ry)" >> "$RESULTS_A"

printf "\n${BOLD}%-14s %s${NC}\n" "a (Å)" "E_total (Ry)" >&2
echo "--------------------------------" >&2

for i in $(seq -$A_NSTEPS $A_NSTEPS); do
    a=$(awk -v c="$A_CENTER" -v s="$A_STEP" -v i="$i" \
        'BEGIN{printf "%.6f", c + i*s}')
    infile="pw_a${a}.in"
    outfile="pw_a${a}.out"

    echo -e "${CYAN}  → a = ${a} Å${NC}" >&2
    write_input "$a" "$C_FIXED" "$infile"

    energy=$(run_pw "$infile" "$outfile")
    if [ $? -ne 0 ] || [ -z "$energy" ]; then
        printf "%-14s %s\n" "$a" "ERROR" >&2
        continue
    fi

    printf "%-14s %s\n" "$a" "$energy" >&2
    echo "$a  $energy" >> "$RESULTS_A"
done

# Encontrar mínimo en a
read min_a min_ea <<< $(find_min "$RESULTS_A")
echo -e "\n  ${GREEN}${BOLD}Mínimo en a = ${min_a} Å  →  E = ${min_ea} Ry${NC}" >&2
A_FIXED=$min_a

# ===========================================================================
# PASO 2: barrido de c (a fijo al mínimo)
# ===========================================================================
echo -e "\n${BOLD}[2/2] Barrido de c  (a fijo = ${A_FIXED} Å)${NC}" >&2
echo "# Barrido de c — a fijo = ${A_FIXED} A   ecutwfc=${ECUTWFC} Ry" > "$RESULTS_C"
echo "# c(Ang)      E_total(Ry)" >> "$RESULTS_C"

printf "\n${BOLD}%-14s %s${NC}\n" "c (Å)" "E_total (Ry)" >&2
echo "--------------------------------" >&2

for i in $(seq -$C_NSTEPS $C_NSTEPS); do
    c=$(awk -v c="$C_CENTER" -v s="$C_STEP" -v i="$i" \
        'BEGIN{printf "%.2f", c + i*s}')
    infile="pw_c${c}.in"
    outfile="pw_c${c}.out"

    echo -e "${CYAN}  → c = ${c} Å${NC}" >&2
    write_input "$A_FIXED" "$c" "$infile"

    energy=$(run_pw "$infile" "$outfile")
    if [ $? -ne 0 ] || [ -z "$energy" ]; then
        printf "%-14s %s\n" "$c" "ERROR" >&2
        continue
    fi

    printf "%-14s %s\n" "$c" "$energy" >&2
    echo "$c  $energy" >> "$RESULTS_C"
done

# Encontrar mínimo en c
read min_c min_ec <<< $(find_min "$RESULTS_C")
echo -e "\n  ${GREEN}${BOLD}Mínimo en c = ${min_c} Å  →  E = ${min_ec} Ry${NC}" >&2

# ===========================================================================
# RESUMEN FINAL
# ===========================================================================
echo "" >&2
echo "================================================================" >&2
echo -e "${BOLD}RESUMEN FINAL${NC}" >&2
echo "----------------------------------------------------------------" >&2
echo -e "  ${GREEN}${BOLD}a_opt = ${min_a} Å${NC}" >&2
echo -e "  ${GREEN}${BOLD}c_opt = ${min_c} Å${NC}" >&2
echo "" >&2
echo "  Usa estos valores en tu SCF / DFPT:" >&2
echo "" >&2
echo "  CELL_PARAMETERS (angstrom)" >&2
awk -v a="$min_a" 'BEGIN{
    printf "    %.10f   0.0000000000   0.0000000000\n", a
    printf "    %.10f   %.10f   0.0000000000\n", -a/2, a*sqrt(3)/2
}' >&2
echo -e "    0.0000000000   0.0000000000   ${min_c}" >&2
echo "" >&2
echo "  Datos barrido a : $RESULTS_A" >&2
echo "  Datos barrido c : $RESULTS_C" >&2
echo "  Log completo    : $LOG_FILE" >&2
echo "================================================================" >&2
