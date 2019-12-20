i=0
set xlabel 'real time'
set xrange [0:10]
do for [i=0:59]{
set output sprintf('pictures/springs/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('data/springs/length%d.dat',i) using 2:3 w l title sprintf('l_{%d}(t)',i)
}
set ylabel 'bias'
set output 'springs/bais.png'
plot 'data/springs/bias.dat' using 2:3 w l notitle
