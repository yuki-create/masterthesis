cd `dirname $0` # change current directry to [ms]
gcc -I/usr/local/include -c ms.c -o ms.o
gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#echo comple completed.
str1="."
M=12
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
mkdir -p ${str1}/results
mkdir -p ${str1}/results/convergence # outputs.dat, parameters and results
#echo simulating ms...
./ms > log_before.txt #./ms N NSTEP
<< COMMENTOUT
gnuplot -persist <<-EOFMarker
set terminal png

set output './results/l0.png'
set xrange [0:11000]
plot './log_before.txt' using 1:2 w p title 'l_0', './log_before.txt' using 1:3 w p title 'ld_0'

# data/outputs.dat NARMAモデルと近似波形のプロット
set xrange [6000:11000]
set ylabel 'y(t)'
set output './results/NARMA2.png'
plot './results/outputs.dat' using 1:2 w l title "ms-NARMA2", './results/outputs.dat' using 1:5 w l title "NARMA2"
set output './results/NARMA10.png'
plot './results/outputs.dat' using 1:3 w l title "ms-NARMA10", './results/outputs.dat' using 1:6 w l title "NARMA10"
set output './results/NARMA20.png'
plot './results/outputs.dat' using 1:4 w l title "ms-NARMA20", './results/outputs.dat' using 1:7 w l title "NARMA20"

set xrange [6000:6500]
set output './results/NARMA2_zoom.png'
plot './results/outputs.dat' using 1:2 w l title "ms-NARMA2", './results/outputs.dat' using 1:5 w l title "NARMA2"
set output './results/NARMA10_zoom.png'
plot './results/outputs.dat' using 1:3 w l title "ms-NARMA10", './results/outputs.dat' using 1:6 w l title "NARMA10"
set output './results/NARMA20_zoom.png'
plot './results/outputs.dat' using 1:4 w l title "ms-NARMA20", './results/outputs.dat' using 1:7 w l title "NARMA20"
unset xrange

# data/springs/length[i].dat のプロット(ずらさない系lとずらした系l_dの比較用)
set xrange [0:0.5]
set yrange [0.9999999:1.0000001]
do for [i=0:$M-1]{
set out sprintf('./results/convergence/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('./springs/length%d.dat',i) using 2:3 w l title sprintf('l_{%d}(t)',i), sprintf('./springs/length%d.dat',i) using 2:4 w l title sprintf('dif-l_{%d}(t) (index:0)',i)
}

# plot lyapunov exponent
set out './results/le.png'
  unset yrange
  set xrange [0:11000] ;
  plot "./results/le.dat" with l
  exit ;
EOFMarker
COMMENTOUT
