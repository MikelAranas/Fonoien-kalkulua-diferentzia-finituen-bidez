#!/usr/bin/env bash
set -e
mkdir -p tmp
NP=12
start=$(date +%s)
mpirun -np $NP pw.x   -in c.scf.in   > c.scf.out
mpirun -np $NP ph.x   -in c.ph.in    > c.ph.out
mpirun -np $NP q2r.x  -in c.q2r.in   > c.q2r.out
mpirun -np $NP matdyn.x -in c.matdyn.in > c.matdyn.out
end=$(date +%s)
echo "Tiempo total: $((end-start)) s"
