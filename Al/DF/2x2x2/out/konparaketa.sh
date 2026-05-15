#!/bin/bash

gnuplot << 'EOF'

set terminal pdf
set output "DF_VS_DFPT_Al.pdf"

set title "DF vs DFPT"
set xlabel "q"
set ylabel "w(q) (cm⁻¹)"

set key top right

# DFPT → línea negra
set style line 1 lc rgb "black" lw 3

# IFC → cruces MUY pequeñas
set style line 11 lc rgb "red" pt 0.2 ps 0.4

plot \
  "al.freq.gp"    using 1:2 with lines ls 1 title "DFPT", \
  "al.freq.gp"    using 1:3 with lines ls 1 notitle, \
  "al.freq.gp"    using 1:4 with lines ls 1 notitle, \
  "Al_DF_bandak.dat" using 1:2 with points ls 11 title "IFC", \
  "Al_DF_bandak.dat" using 1:3 with points ls 11 notitle, \
  "Al_DF_bandak.dat" using 1:4 with points ls 11 notitle

EOF

echo "Plot guardado en DF_VS_DFPT_Al.pdf"
