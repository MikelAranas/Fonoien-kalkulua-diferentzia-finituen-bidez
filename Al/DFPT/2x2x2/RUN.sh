#!/usr/bin/env bash
set -e

NP=12

start=$(date +%s)

mpirun -np $NP pw.x     -in al.scf.in    > al.scf.out
mpirun -np $NP ph.x     -in al.ph.in     > al.ph.out
mpirun -np $NP q2r.x    -in al.q2r.in    > al.q2r.out
mpirun -np $NP matdyn.x -in al.matdyn.in > al.matdyn.out

end=$(date +%s)
echo "Tiempo total: $((end-start)) s"
