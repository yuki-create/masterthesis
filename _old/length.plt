i=0
set terminal png
set xlabel 'real time'
do for [i=0:$M-1]{
set output sprintf('$dirname2/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('$dirname1/springs/length%d.dat',i) using 2:3 w l title sprintf('l_{%d}(t)',i)
}
set terminal aqua
set output
