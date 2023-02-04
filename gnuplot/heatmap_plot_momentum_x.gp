#! /opt/homebrew/bin/gnuplot

cd 'data'

set terminal gif animate delay 5
set output 'gif/momentumX_Plot.gif'
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "momentumX(x,y)"
stats 'test_momentumX_Results.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot 'test_momentumX_Results.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}