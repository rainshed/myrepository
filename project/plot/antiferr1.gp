unset key
set xlabel "X"
set ylabel "Y"
plot "data/antiferr1.dat" with image
set term pngcairo
set output "picture/antiferr1.png"
replot
set output

