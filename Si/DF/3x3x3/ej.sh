#!/bin/bash
set -euo pipefail

# =============================================================================
# ej.sh — versión PARALELA POR LOTES del lanzador de desplazados
#
# Igual que el ej.sh secuencial (mismo run_one, mismo pie TIMING en cada .out,
# mismo scratch aislado por job), pero lanza PARALLEL jobs a la vez y espera
# al lote completo antes del siguiente. Con PARALLEL=4 y NP=12 → 48 cores.
#
# La métrica de tiempos NO cambia: denborak.py lee el WALL del PWSCF de cada
# .out, que es por-job. Solo se reduce el tiempo de pared total (~/PARALLEL).
#
# ⚠ MEMORIA: en 6x6x6 (432 átomos) vigilar la RAM con 4 jobs simultáneos.
#   Si se queda corto: PARALLEL=2 NP=24 bash ej.sh
# ⚠ ANCLAJE DE PROCESOS (pinning): con 4 mpirun simultáneos hay que evitar
#   que cada MPI ancle sus 12 procesos a los MISMOS cores. Cómo se desactiva
#   depende del MPI:
#     - Intel MPI (mpiexec Hydra, el del servidor): export I_MPI_PIN=0
#       (ya se hace abajo; mpiexec de Intel NO entiende --bind-to)
#     - OpenMPI: lanzar con  MPIRUN_FLAGS="--bind-to none" bash ej.sh
# =============================================================================

# --- Config ---
NP=${NP:-12}                          # cores MPI por job
PARALLEL=${PARALLEL:-4}               # jobs simultáneos (NP*PARALLEL <= cores)
MPIRUN_FLAGS=${MPIRUN_FLAGS:-}        # vacío por defecto (Intel MPI no traga --bind-to)
INPUT_DIR="in"
OUTPUT_DIR="out"
SCRATCH_ROOT="./scratch_runs"
BASE_SCF_OUTDIR="./tmp"               # contiene .save del SCF base
PREFIX="si_sup"
PW_COMMAND="pw.x"
LOG_DIR="./logs"

mkdir -p "$OUTPUT_DIR" "$SCRATCH_ROOT" "$LOG_DIR"
export OMP_NUM_THREADS=1
export I_MPI_PIN=${I_MPI_PIN:-0}      # Intel MPI: sin anclaje, el SO reparte cores
ulimit -s unlimited

# --- Utilidades de tiempo ---
fmt_duration() {
  local s=$1
  printf "%02d:%02d:%02d" $((s/3600)) $(((s%3600)/60)) $((s%60))
}

run_one() {
  local infile="$1"
  local base="$(basename "$infile" .in)"
  local scratch="$SCRATCH_ROOT/$base"
  local out="$OUTPUT_DIR/${base}.out"
  local log="$LOG_DIR/${base}.log"
  local tmpin="$scratch/${base}.in"

  local t0=$(date +%s)
  local t0_iso=$(date -Iseconds)

  echo "[PREP] $base (NP=$NP) @ $t0_iso"
  rm -rf "$scratch"; mkdir -p "$scratch"

  # Copia limpia del save base al scratch del job
  rsync -a --delete "$BASE_SCF_OUTDIR/$PREFIX.save/" "$scratch/$PREFIX.save/" >>"$log" 2>&1

  # Inyecta el outdir de scratch en el input temporal
  sed -e "s|outdir *= *'.*'|outdir = '${scratch}'|" "$infile" > "$tmpin"

  echo "[RUN ] $base"
  if mpirun $MPIRUN_FLAGS -np "$NP" "$PW_COMMAND" -in "$tmpin" > "$out" 2>>"$log"; then
    :
  else
    echo "[FAIL] $base — revisa $log"
  fi

  local t1=$(date +%s)
  local t1_iso=$(date -Iseconds)
  local dt=$((t1 - t0))

  {
    echo ""
    echo "===== TIMING ====="
    echo "Job:           $base"
    echo "Start:         $t0_iso"
    echo "End:           $t1_iso"
    echo "Duration (s):  $dt"
    echo "Duration (H:M:S): $(fmt_duration "$dt")"
    echo "=================="
  } >> "$out"

  echo "[DONE] $base @ $t1_iso  (duración: $(fmt_duration "$dt"))" >> "$log"
  echo "[DONE] $base  (duración: $(fmt_duration "$dt"))"
}

# --- Chequeos y lista de inputs ---
shopt -s nullglob
inputs=( "$INPUT_DIR"/*.in )
shopt -u nullglob
if [ ${#inputs[@]} -eq 0 ]; then
  echo "ERROR: No hay .in en $INPUT_DIR"; exit 1
fi
[ -d "$BASE_SCF_OUTDIR/$PREFIX.save" ] || {
  echo "ERROR: No existe $BASE_SCF_OUTDIR/$PREFIX.save/"; exit 1; }

# --- Cronómetro global ---
global_t0=$(date +%s)
global_t0_iso=$(date -Iseconds)
total_jobs=${#inputs[@]}
echo "Lanzando $total_jobs desplazados (lotes de $PARALLEL, NP=$NP por job) @ $global_t0_iso"

# --- Bucle por lotes ---
job_idx=0
batch_idx=0
i=0
while [ $i -lt $total_jobs ]; do
  batch_idx=$((batch_idx + 1))
  batch=( "${inputs[@]:$i:$PARALLEL}" )
  echo ""
  echo ">>> [LOTE $batch_idx] ${#batch[@]} jobs:" \
       "$(for f in "${batch[@]}"; do basename "$f" .in; done | tr '\n' ' ')"

  pids=()
  for infile in "${batch[@]}"; do
    run_one "$infile" &
    pids+=( $! )
  done

  # Espera del lote completo; si algún job revienta, se informa y se sigue
  fail=0
  for pid in "${pids[@]}"; do
    wait "$pid" || fail=1
  done
  [ $fail -ne 0 ] && echo "[AVISO] Algún job del lote $batch_idx falló — revisa $LOG_DIR"

  i=$((i + PARALLEL))
  job_idx=$((i < total_jobs ? i : total_jobs))
  now=$(date +%s)
  elapsed=$((now - global_t0))
  echo "    Progreso: $job_idx/$total_jobs — tiempo transcurrido: $(fmt_duration "$elapsed")"
done

global_t1=$(date +%s)
global_dt=$((global_t1 - global_t0))

echo ""
echo "Todos los desplazados terminados."
echo "Inicio:     $global_t0_iso"
echo "Fin:        $(date -Iseconds)"
echo "Total (pared): $(fmt_duration "$global_dt")"
echo "Outputs:    $OUTPUT_DIR"
echo "Logs:       $LOG_DIR"
echo ""
echo "Recordatorio: t_DF del informe = suma de WALL de cada .out (denborak.py),"
echo "no este tiempo de pared."

