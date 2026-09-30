set terminal pngcairo size 800,800 enhanced font 'Helvetica,10'
set output 'pic.png'
set datafile separator ","
set autoscale fix

plot "plot_data.csv" using 1:2:3 with points pt 7 ps 0.5 lc rgb variable notitle
