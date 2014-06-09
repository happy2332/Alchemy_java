set term postscript eps enhanced color 'Helvetica' 20
#set title "Number of Launched and Visited Mobiles"
set output "data.eps"
#Launched = "#99ffff"; Visited = "#4671d5";
myColor = "#99ffff"
set auto x
set auto y
set yrange [0:25000]
set y2range [0:150]
#set xtics nomirror rotate by -45 font ",12"
set xtics nomirror font ",16"
set ytics font ",16"
set y2tics font ",16"
#set style data histogram
#set style histogram cluster gap 1
#set style fill solid border -1
#set boxwidth 0.9
set style line 1 lt 1 lw 3 pt 3 linecolor rgb "red"
set style line 2 lt 1 lw 3 pt 3 linecolor rgb "blue"
set xtic scale 1
     set xdata time
     set timefmt "%d %b"
     set format x "%d %b"
     set xrange ["8 jun":"4 july"]
     set xtics ("8 jun", "12 jun", "30 jun","4 july")


set xlabel "Date" font "12" offset 0,0
set ylabel "Price" font "12" offset 2,0
set y2label "No. of converts" font "10" offset -2,0
set datafile separator "\t"
# 2, 3, 4, 5 are the indexes of the columns; 'fc' stands for 'fillcolor'
#plot 'convert_month.data' using 2:xtic(1) ti col fc rgb Launched, '' u 3 ti col fc rgb Visited
plot 'sony_xperia_price_convert.data' using 1:2 ls 1 w l ti "Price",\
'sony_xperia_price_convert.data' using 1:3 ls 2 w l ti "No. of Converts" axes x1y2
