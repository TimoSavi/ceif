#!/bin/sh

cat plot_data.csv | awk "-F," '$4 > 0.499 && $4 < 0.501' > contour_d50.csv
cat plot_data.csv | awk "-F," '$4 > 0.549 && $4 < 0.551' > contour_d55.csv
cat plot_data.csv | awk "-F," '$4 > 0.599 && $4 < 0.601' > contour_d60.csv
cat plot_data.csv | awk "-F," '$4 > 0.699 && $4 < 0.701' > contour_d70.csv
cat plot_data.csv | awk "-F," '$4 > 0.799 && $4 < 0.801' > contour_d80.csv

