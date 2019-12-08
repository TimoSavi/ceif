set terminal png
set output 'pic.png'
set xrange [0:40]
set yrange [0:40]
set datafile separator ","

plot "plot_data.csv" using 1:2:3 with points pt 7 ps 0.5 lc rgb variable notitle, \
      datafile using 1:2 with points pt 7 ps 0.25 lc rgb "black" notitle, \
     "contour_d50.csv" using 1:2 with points pt 7 ps 0.25 lc rgb "blue" notitle, \
     "contour_d55.csv" using 1:2 with points pt 7 ps 0.25 lc rgb "brown" notitle, \
     "contour_d60.csv" using 1:2 with points pt 7 ps 0.25 lc rgb "green" notitle, \
     "contour_d70.csv" using 1:2 with points pt 7 ps 0.25 lc rgb "yellow" notitle, \
     "contour_d80.csv" using 1:2 with points pt 7 ps 0.25 lc rgb "red" notitle

