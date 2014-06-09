set term postscript eps enhanced color 'Helvetica' 20
set output "data.eps"
set logscale y
set xlabel "Time in seconds"
set ylabel "False weight"
plot "data1.data" using 1:2:3:4 with yerrorlines title "data1",\
 "data2.data" using 1:2:3:4 with yerrorlines title "data2"
