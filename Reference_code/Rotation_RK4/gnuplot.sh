gnuplot <<- EOF
	reset;
	set terminal postscript eps color font "Times-Roman,26" enhanced;
	set output "Graph.eps";
set grid;
	set key center top outside horizontal
	set xlabel "time [s]"
	set ylabel "angular velocity [rad/s]"
	plot "output.dat" using 1:2 with lines dt 1 lw 4 lt rgb "black" title "{/Symbol w}_1","output.dat" using 1:3 with lines dt 1 lw 4 lt rgb "red" title "{/Symbol w}_2","output.dat" using 1:4 with lines dt 1 lw 4 lt rgb "green" title "{/Symbol w}_3";
	set ou;

EOF

#convert -density 500 Graph.eps Graph.jpg
#rm Graph.eps

