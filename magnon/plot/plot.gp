set xlabel "X"
set ylabel "E"

plot "data/energyband.dat" u 1:2,"data/energyband.dat" u 1:3,\
	     "data/energyband.dat" u 1:4,"data/energyband.dat" u 1:5,"data/energyband.dat" u 1:6,"data/energyband.dat" u 1:7,\
# "data/energyband.dat" u 1:8,"data/energyband.dat" u 1:9,"data/energyband.dat" u 1:10,"data/energyband.dat" u 1:11,\
#	     "data/energyband.dat" u 1:12,"data/energyband.dat" u 1:13

set term pngcairo
set output "plot/energyband.png"
replot
set output
