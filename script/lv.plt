# set term postscript eps size 1280,720
set term svg size 1280,720
set out "lv.svg"

set title "Original function"

set autoscale

set xrange [ 0 : 200 ]
set autoscale y

set xlabel "Time"
set ylabel "Count"

plot "lv.dat" using 1:2 w lines, \
     "lv.dat" using 1:3 w lines
