# Plotting script for complex2d dataset and CEIF heatmap
# Usage:
#   1. Train and generate score map:
#      ./src/ceif -l test/complex2d.csv -T 0.5 -i 300 -O 0s -p "%d,0x%x" -o test/plot_complex2d.csv
#   2. Render plots:
#      gnuplot test/plot_complex2d.gp

set datafile separator ','

# 1. Render heatmap
set terminal pngcairo size 800,800 enhanced font 'Helvetica,10'
set output 'test/complex2d.png'
set title 'CEIF Anomaly Score Heatmap - Complex 2D Test Dataset'
set xlabel 'Dimension 1 (Scale: 1000 - 9000)'
set ylabel 'Dimension 2 (Scale: 0.2 - 1.8)'
plot 'test/plot_complex2d.csv' using 1:2:3 with points pt 7 ps 0.4 lc rgb variable notitle

# 2. Render raw points
set output 'test/complex2d_points.png'
set title 'Raw 2D Synthetic Test Data (1465 points)'
plot 'test/complex2d.csv' using 1:2 with points pt 7 ps 0.6 lc rgb '#0055aa' notitle
