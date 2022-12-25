set terminal jpeg
set output 'pic.jpg'
set datafile separator ","

plot "plot_data.csv" using 1:2:3 with points pt 7 ps 0.5 lc rgb variable notitle
