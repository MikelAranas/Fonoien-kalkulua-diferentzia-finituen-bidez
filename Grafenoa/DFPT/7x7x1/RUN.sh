#!/usr/bin/env bash
# =============================================================================
# RUN.sh — Dispersión de fonones del grafeno
# pw.x → ph.x → q2r.x → matdyn.x → bandak_cm.py
# =============================================================================

set -euo pipefail

NP=12
TOTAL_STEPS=5

# ---------------------------------------------------------------------------
# Colores
# ---------------------------------------------------------------------------
RED='\033[0;31m'
GREEN='\033[0;32m'
CYAN='\033[0;36m'
YELLOW='\033[1;33m'
BOLD='\033[1m'
NC='\033[0m'

# ---------------------------------------------------------------------------
# Funciones
# ---------------------------------------------------------------------------
step_start() {
    local n=$1
    local name=$2
    echo ""
    echo -e "${CYAN}${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "${CYAN}${BOLD}  Paso ${n}/${TOTAL_STEPS}: ${name}${NC}"
    echo -e "${CYAN}${BOLD}━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━━${NC}"
    echo -e "  Inicio : $(date '+%H:%M:%S')"
    STEP_START=$(date +%s)
}

step_end() {
    local n=$1
    local name=$2
    local outfile=$3
    local elapsed=$(( $(date +%s) - STEP_START ))
    echo -e "  Fin    : $(date '+%H:%M:%S')  (${elapsed} s)"
    # Muestra la última línea relevante del output como confirmación
    local last
    last=$(grep -E "JOB DONE|Total|total energy|Diagonalizing|q = " "$outfile" 2>/dev/null | tail -1 || true)
    if [ -n "$last" ]; then
        echo -e "  ${GREEN}✓ ${last// /  }${NC}"
    fi
    echo -e "  ${GREEN}${BOLD}Paso ${n} completado.${NC}"
}

check_input() {
    local f=$1
    if [ ! -f "$f" ]; then
        echo -e "${RED}ERROR: no se encuentra el fichero de input '${f}'. Abortando.${NC}"
        exit 1
    fi
}

# ---------------------------------------------------------------------------
# Cabecera
# ---------------------------------------------------------------------------
echo -e "${BOLD}"
echo "============================================================"
echo "   Grafeno — Dispersión de fonones con QE"
echo "   $(date '+%Y-%m-%d %H:%M:%S')"
echo "   Cores MPI: ${NP}"
echo "============================================================"
echo -e "${NC}"

# Verificar inputs antes de empezar
echo -e "${YELLOW}Verificando ficheros de input...${NC}"
for f in c.scf.in c.ph.in c.q2r.in c.matdyn.in; do
    check_input "$f"
    echo -e "  ${GREEN}✓${NC} ${f}"
done
echo ""

mkdir -p tmp
GLOBAL_START=$(date +%s)

# ===========================================================================
# PASO 1: SCF
# ===========================================================================
step_start 1 "SCF (pw.x)"
mpirun -np $NP pw.x -in c.scf.in > c.scf.out
step_end 1 "SCF" c.scf.out

# Extraer energía y parámetro de red del SCF
ETOT=$(grep "!    total energy" c.scf.out | tail -1 | awk '{print $5, $6}')
echo -e "  Energía total : ${ETOT}"

# ===========================================================================
# PASO 2: DFPT — fonones (ph.x)
# ===========================================================================
step_start 2 "DFPT fonones (ph.x)"
echo -e "  ${YELLOW}Este paso puede tardar varios minutos...${NC}"
mpirun -np $NP ph.x -in c.ph.in > c.ph.out
step_end 2 "DFPT" c.ph.out

# Número de q-points calculados
NQ=$(grep "Dynamical matrix file" c.ph.out 2>/dev/null | wc -l || echo "?")
echo -e "  Matrices dinámicas escritas: ${NQ}"

# ===========================================================================
# PASO 3: Transformada de Fourier (q2r.x)
# ===========================================================================
step_start 3 "Fourier inversa — constantes de fuerza (q2r.x)"
mpirun -np $NP q2r.x -in c.q2r.in > c.q2r.out
step_end 3 "q2r" c.q2r.out

if [ -f graphene.fc ]; then
    echo -e "  ${GREEN}✓ graphene.fc generado${NC}  ($(du -h graphene.fc | cut -f1))"
else
    echo -e "  ${RED}AVISO: no se encontró graphene.fc${NC}"
fi

# ===========================================================================
# PASO 4: Dispersión de fonones (matdyn.x)
# ===========================================================================
step_start 4 "Dispersión de fonones (matdyn.x)"
mpirun -np $NP matdyn.x -in c.matdyn.in > c.matdyn.out
step_end 4 "matdyn" c.matdyn.out

if [ -f graphene.freq ]; then
    echo -e "  ${GREEN}✓ graphene.freq generado${NC}"
    # Avisar si hay frecuencias negativas
    NEG=$(grep -c "^\s*-" graphene.freq.gp 2>/dev/null || echo 0)
    if [ "$NEG" -gt 0 ]; then
        echo -e "  ${YELLOW}⚠ Hay ${NEG} frecuencias negativas en graphene.freq.gp — revisar${NC}"
    else
        echo -e "  ${GREEN}✓ Sin frecuencias negativas detectadas${NC}"
    fi
fi

# ===========================================================================
# PASO 5: Plot (bandak_cm.py)
# ===========================================================================
step_start 5 "Generando figura (bandak_cm.py)"
python3 bandak_cm.py \
    c.matdyn.out c.ph.out c.q2r.out c.scf.out \
    graphene_bandak_cm.dat \
    graphene.dyn0 graphene.dyn2 graphene.dyn4 \
    graphene.freq matdyn.modes path.dat \
    > bandak_cm.log 2>&1
step_end 5 "bandak_cm.py" bandak_cm.log

if [ -f graphene_bandak.pdf ]; then
    echo -e "  ${GREEN}${BOLD}✓ Figura generada: graphene_bandak.pdf${NC}"
fi

# ===========================================================================
# RESUMEN FINAL
# ===========================================================================
TOTAL=$(( $(date +%s) - GLOBAL_START ))
MINS=$(( TOTAL / 60 ))
SECS=$(( TOTAL % 60 ))

echo ""
echo -e "${GREEN}${BOLD}"
echo "============================================================"
echo "   ✓ CÁLCULO COMPLETADO"
printf "   Tiempo total: %dm %02ds\n" $MINS $SECS
echo "============================================================"
echo -e "${NC}"
echo "  Ficheros generados:"
for f in c.scf.out c.ph.out c.q2r.out c.matdyn.out \
          graphene.fc graphene.freq graphene_bandak.pdf; do
    if [ -f "$f" ]; then
        echo -e "    ${GREEN}✓${NC} ${f}  ($(du -h "$f" | cut -f1))"
    else
        echo -e "    ${RED}✗${NC} ${f}  — no encontrado"
    fi
done
echo ""
