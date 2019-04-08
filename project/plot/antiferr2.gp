unset key
set xlabel "X"
set ylabel "Y"
plot "data/antiferr2.dat" with image
set term pngcairo
set output "picture/antiferr2.png"
replot
set output

