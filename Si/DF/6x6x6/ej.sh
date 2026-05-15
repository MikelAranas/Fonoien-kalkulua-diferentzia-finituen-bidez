#!/bin/bash
set -euo pipefail

# --- Config ---
NP=${NP:-12}                          # cores MPI por job (por defecto 12)
INPUT_DIR="in"          # carpeta de .in desplazados
OUTPUT_DIR="out"            # carpeta de .out finales
SCRATCH_ROOT="./scratch_runs"         # scratch aislado por job
BASE_SCF_OUTDIR="./tmp" # contiene .save del SCF base
PREFIX="si_sup"
PW_COMMAND="pw.x"
LOG_DIR="./logs"

mkdir -p "$OUTPUT_DIR" "$SCRATCH_ROOT" "$LOG_DIR"
export OMP_NUM_THREADS=1
ulimit -s unlimited

# --- Utilidades de tiempo ---
fmt_duration() {  # formatea segundos a H:M:S
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

  # Marcas temporales (inicio job)
  local t0=$(date +%s)
  local t0_iso=$(date -Iseconds)

  echo "[PREP] $base (NP=$NP) @ $t0_iso"
  rm -rf "$scratch"; mkdir -p "$scratch"

  # Copia limpia del save base al scratch del job
  rsync -a --delete "$BASE_SCF_OUTDIR/$PREFIX.save/" "$scratch/$PREFIX.save/" >>"$log" 2>&1

  # Inyecta el outdir de scratch en el input temporal
  sed -e "s|outdir *= *'.*'|outdir = '${scratch}'|" "$infile" > "$tmpin"

  echo "[RUN ] $base"
  # Ejecuta QE; stdout al .out, stderr al .log
  if mpirun -np "$NP" "$PW_COMMAND" -in "$tmpin" > "$out" 2>>"$log"; then
    :
  else
    echo "[FAIL] $base — revisa $log"
  fi

  # Marcas temporales (fin job)
  local t1=$(date +%s)
  local t1_iso=$(date -Iseconds)
  local dt=$((t1 - t0))

  # Añade una cabecera de timing al out y al log
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

  {
    echo "[DONE] $base @ $t1_iso  (duración: $(fmt_duration "$dt"))"
  } >> "$log"

  echo "[DONE] $base  (duración: $(fmt_duration "$dt"))"
}

# --- Chequeos y preparación de lista de inputs ---
shopt -s nullglob
inputs=( "$INPUT_DIR"/*.in )
shopt -u nullglob
if [ ${#inputs[@]} -eq 0 ]; then
  echo "ERROR: No hay .in en $INPUT_DIR"; exit 1
fi
[ -d "$BASE_SCF_OUTDIR/$PREFIX.save" ] || {
  echo "ERROR: No existe $BASE_SCF_OUTDIR/$PREFIX.save/"; exit 1; }

# Cronómetro global
global_t0=$(date +%s)
global_t0_iso=$(date -Iseconds)
echo "Lanzando cálculos desplazados (secuencial, NP=$NP) @ $global_t0_iso"

# Bucle secuencial
job_idx=0
total_jobs=${#inputs[@]}
for infile in "${inputs[@]}"; do
  job_idx=$((job_idx + 1))
  echo ">>> [JOB $job_idx/$total_jobs] $(basename "$infile" .in)"
  run_one "$infile"

  # Informe de progreso tras cada job
  now=$(date +%s)
  elapsed=$((now - global_t0))
  echo "    Progreso: $job_idx/$total_jobs — tiempo transcurrido: $(fmt_duration "$elapsed")"
done

global_t1=$(date +%s)
global_dt=$((global_t1 - global_t0))
global_t1_iso=$(date -Iseconds)

echo ""
echo "Todos los desplazados terminados."
echo "Inicio:     $global_t0_iso"
echo "Fin:        $global_t1_iso"
echo "Total:      $(fmt_duration "$global_dt")"
echo "Outputs:    $OUTPUT_DIR"
echo "Logs:       $LOG_DIR"


