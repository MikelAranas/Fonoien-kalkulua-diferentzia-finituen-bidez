#!/bin/bash

gnuplot << 'EOF'

set terminal pdf
set output "DF_VS_DFPT_Grafenoa.pdf"

set title "DF vs DFPT (Grafenoa)"
set xlabel "q"
set ylabel "w(q) (cm^{-1})"

set key top right

# DFPT → línea negra
set style line 1 lc rgb "black" lw 3

# IFC → cruces MUY pequeñas
set style line 11 lc rgb "red" pt 0.2 ps 0.4

plot \
  "graphene.freq.gp" using 1:2 with lines  ls 1 title "DFPT", \
  "graphene.freq.gp" using 1:3 with lines  ls 1 notitle, \
  "graphene.freq.gp" using 1:4 with lines  ls 1 notitle, \
  "graphene.freq.gp" using 1:5 with lines  ls 1 notitle, \
  "graphene.freq.gp" using 1:6 with lines  ls 1 notitle, \
  "graphene.freq.gp" using 1:7 with lines  ls 1 notitle, \
  "Grafenoa_DF_bandak.dat" using 1:2 with points ls 11 title "IFC", \
  "Grafenoa_DF_bandak.dat" using 1:3 with points ls 11 notitle, \
  "Grafenoa_DF_bandak.dat" using 1:4 with points ls 11 notitle, \
  "Grafenoa_DF_bandak.dat" using 1:5 with points ls 11 notitle, \
  "Grafenoa_DF_bandak.dat" using 1:6 with points ls 11 notitle, \
  "Grafenoa_DF_bandak.dat" using 1:7 with points ls 11 notitle

EOF

echo "Plot guardado en DF_VS_DFPT_Grafenoa.pdf"

