#! /opt/homebrew/bin/gnuplot

cd 'data'

set terminal gif animate delay 5
set output 'gif/rhoPlot.gif'
set hidden3d
stats 'test_rhoResults.dat' using 4 nooutput
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "rho(x,y)"
set zrange [STATS_min:STATS_max]
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    splot 'test_rhoResults.dat' index (i-1) using 3:2:4 with lines palette notitle
}

set terminal gif animate delay 5
set output 'gif/momentumX_Plot.gif'
set hidden3d
stats 'test_momentumX_Results.dat' using 4 nooutput
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "momentumX(x,y)"
set zrange [STATS_min:STATS_max]
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    splot 'test_momentumX_Results.dat' index (i-1) using 2:3:4 with lines palette notitle
}

set terminal gif animate delay 5
set output 'gif/momentumY_Plot.gif'
set hidden3d
stats 'test_momentumY_Results.dat' using 4 nooutput
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "momentumY(x,y)"
set zrange [STATS_min:STATS_max]
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    splot 'test_momentumY_Results.dat' index (i-1) using 2:3:4 with lines palette notitle
}

set terminal gif animate delay 5
set output 'gif/energyPlot.gif'
set hidden3d
stats 'test_energyResults.dat' using 4 nooutput
set xyplane at 0
set xlabel "x"
set ylabel "y"
set zlabel "energy(x,y)"
set zrange [STATS_min:STATS_max]
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    splot 'test_energyResults.dat' index (i-1) using 2:3:4 with lines palette notitle
}
