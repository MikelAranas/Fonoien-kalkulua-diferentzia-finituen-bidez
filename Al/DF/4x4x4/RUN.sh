#!/usr/bin/env bash
set -e
NP=12
start=$(date +%s)

mpirun -np $NP pw.x -in al_sup.scf.in > al_sup.scf.out
python crearinputs.py
bash ej.sh
(cd out && python IFC_murriztue.py)
(cd out && python D.py)

end=$(date +%s)
echo "Tiempo total: $((end-start)) s"

