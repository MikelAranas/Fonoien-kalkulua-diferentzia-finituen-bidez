#!/usr/bin/env bash
set -e

NP=12

start=$(date +%s)

mpirun -np $NP pw.x     -in si.scf.in    > si.scf.out
mpirun -np $NP ph.x     -in si.ph.in     > si.ph.out
mpirun -np $NP q2r.x    -in si.q2r.in    > si.q2r.out
mpirun -np $NP matdyn.x -in si.matdyn.in > si.matdyn.out

end=$(date +%s)
echo "Tiempo total: $((end-start)) s"
