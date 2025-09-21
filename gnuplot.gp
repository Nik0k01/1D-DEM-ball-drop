set output "Plots/Energy_plot.png"
set terminal pngcairo size 800,600 enhanced font "Arial,14"
set grid
set style line 1 lt 1 lw 2 lc rgb "blue"
set style line 2 lt 2 lw 2 lc rgb "red"
set style line 3 lt 3 lw 2 lc rgb "green"
set style line 4 lt 3 lw 2 lc rgb "black"
set style line 5 lt 3 lw 2 lc rgb "pink"
set xlabel "Time (s)"
set ylabel "Energy (J)"
set title "Energy vs Time"
plot "Energy_data.dat" using 1:2 with lines ls 1 title "Potential", "Energy_data.dat" using 1:3 with lines ls 2 title "Kinetic", "Energy_data.dat" using 1:4 with lines ls 3 title "Elastic", "Energy_data.dat" using 1:5 with lines ls 4 title "Dispersed", "Energy_data.dat" using 1:6 with lines ls 5 title "Total",
