M=12
dirname = 'waves'
i=0
set terminal png
set xlabel 'real time'

do for [i=0:M-1]{
set output './'.dirname.sprintf('/l_%02d.png',i)
set ylabel sprintf('l_%d(t)',i)
plot sprintf('./springs/length%d.dat',i) using 2:3 w l title sprintf('l_%d(t)',i)
}
set terminal aqua
set output
