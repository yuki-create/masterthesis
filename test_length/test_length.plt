M=12
n=0
set xlabel 'real time'
set ylabel 'l_0(t)'

do for []
plot './springs/length0.dat' using 2:3 w l title'l_0(t)'
