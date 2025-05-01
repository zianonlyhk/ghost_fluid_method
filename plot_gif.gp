#! /usr/bin/gnuplot

cd 'output'

# Ensure gif directory exists
system("mkdir -p gif")

# Get simulation name from config.txt
sim_name = system("awk '/^runName/ {print $2}' ../config.txt")

# Common plot settings variables
TERM_SETTINGS = "gif animate delay 3"
PALETTE_SETTINGS = "gray"
AXIS_LABELS = 'set xlabel "x"; set ylabel "y"'

# Plot density
set terminal @TERM_SETTINGS
set palette @PALETTE_SETTINGS
set xyplane at 0
eval(AXIS_LABELS)
set output 'gif/'.sim_name.'_rho.gif'
stats sim_name.'_rhoResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot sim_name.'_rhoResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}

# Plot velocity magnitude
set terminal @TERM_SETTINGS
set palette @PALETTE_SETTINGS
set xyplane at 0
eval(AXIS_LABELS)
set output 'gif/'.sim_name.'_velMag.gif'
stats sim_name.'_velMagResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot sim_name.'_velMagResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}

# Plot pressure
set terminal @TERM_SETTINGS
set palette @PALETTE_SETTINGS
set xyplane at 0
eval(AXIS_LABELS)
set output 'gif/'.sim_name.'_pressure.gif'
stats sim_name.'_pressureResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot sim_name.'_pressureResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}

# Plot mockschlieren
set terminal @TERM_SETTINGS
set palette @PALETTE_SETTINGS
set xyplane at 0
eval(AXIS_LABELS)
set output 'gif/'.sim_name.'_mockschlieren.gif'
stats sim_name.'_msResults.dat' using 4 nooutput
set cbrange [STATS_min:STATS_max]
do for [i=1:int(STATS_blocks)-1] {
    plot sim_name.'_msResults.dat' index (i-1) using 2:3:4 with image palette notitle
    set xrange [GPVAL_DATA_X_MIN:GPVAL_DATA_X_MAX]
    set yrange [GPVAL_DATA_Y_MIN:GPVAL_DATA_Y_MAX]
}
