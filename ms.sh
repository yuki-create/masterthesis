cd `dirname $0` # change current directry to [ms]
gcc -I/usr/local/include -c ms.c -o ms.o
gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#echo comple completed.
# [0]./ms, [1]debug_flag [2]N(pow of integer), [3]mu_k, [4]sigma_k, [5]mu_g, [6]sigma_g, [7] dirname
# 正規分布の muは平均、sigmaは標準偏差
time=`date "+%m%d_%H%M%S"`
dirname1="./${time}"
#dirname1="delta-test"
N=6
debug=1
mu_k=1000
sigma_k=300
mu_g=0.01
sigma_g=0.003
d_idx=0
L=$(($N**2)) #点の個数
M=$(($(($(($N-1))*$N))*2)) # ばねの本数
if [ $debug -eq 1 ] ; then
mkdir -p ${dirname1}
mkdir -p ${dirname1}/data #outputs, le .dat
mkdir -p ${dirname1}/data/points
mkdir -p ${dirname1}/data/springs
mkdir -p ${dirname1}/results # parameters and results
mkdir -p ${dirname1}/results/pictures
mkdir -p ${dirname1}/results/pictures/springs
mkdir -p ${dirname1}/results/pictures/convergence
fi
./ms $debug $N $mu_k $sigma_k $mu_g $sigma_g ${dirname1}
#./ms 1 8 1000 300 0.05 0.01 "./test"
dirname2="${dirname1}/results/pictures"
if [ $debug -eq 1 ] ; then

gnuplot -persist <<-EOFMarker
set terminal png
# data/springs/length[i].dat のプロット
i=0
set xlabel 'real time'
set xrange [15:27.5]
do for [i=0:$M-1]{
set output sprintf('$dirname2/springs/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('$dirname1/data/springs/length%d.dat',i) using 2:3 w l title sprintf('l_{%d}(t)',i)
}
set ylabel 'bias'
set output '$dirname2/springs/bais.png'
plot '$dirname1/data/springs/bias.dat' using 2:3 w l notitle


# data/springs/length[i].dat のプロット(ずらさない系lとずらした系l_dの比較用)
set xrange [0:2.5]
do for [i=0:$M-1]{
set out sprintf('$dirname2/convergence/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('$dirname1/data/springs/length%d.dat',i) using 2:3 w l title sprintf('l_{%d}(t)',i), sprintf('$dirname1/data/springs/length%d.dat',i) using 2:4 w l title sprintf('dif-l_{%d}(t) (index:$d_idx)',i)
}

# data/input.dat 入力時系列のプロット
set ylabel 'I(t)'
set xrange [15:16.5]
set output '$dirname2/input.png'
plot '$dirname1/data/input.dat' using 2:3 w l notitle
unset xrange

# data/outputs.dat NARMAモデルと近似波形のプロット
set xrange [15:27.5]
set ylabel 'y(t)'
set output '$dirname2/NARMA2.png'
plot '$dirname1/data/outputs.dat' using 2:3 w l title "ms-NARMA2", '$dirname1/data/outputs.dat' using 2:7 w l title "NARMA2"
set output '$dirname2/NARMA10.png'
plot '$dirname1/data/outputs.dat' using 2:4 w l title "ms-NARMA10", '$dirname1/data/outputs.dat' using 2:8 w l title "NARMA10"
set output '$dirname2/NARMA20.png'
plot '$dirname1/data/outputs.dat' using 2:5 w l title "ms-NARMA20", '$dirname1/data/outputs.dat' using 2:9 w l title "NARMA20"
set output '$dirname2/NARMA30.png'
plot '$dirname1/data/outputs.dat' using 2:6 w l title "ms-NARMA30", '$dirname1/data/outputs.dat' using 2:10 w l title "NARMA30"

set xrange [15:16.5]
set output '$dirname2/NARMA2_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:3 w l title "ms-NARMA2", '$dirname1/data/outputs.dat' using 2:7 w l title "NARMA2"
set output '$dirname2/NARMA10_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:4 w l title "ms-NARMA10", '$dirname1/data/outputs.dat' using 2:8 w l title "NARMA10"
set output '$dirname2/NARMA20_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:5 w l title "ms-NARMA20", '$dirname1/data/outputs.dat' using 2:9 w l title "NARMA20"
set output '$dirname2/NARMA30_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:6 w l title "ms-NARMA30", '$dirname1/data/outputs.dat' using 2:10 w l title "NARMA30"
unset xrange

# results/le.dat　のプロット
set ylabel 'lyapunov exponent'
set output '$dirname2/le.png'
plot '$dirname1/data/le.dat' using 2:3 w l notitle

#
set terminal aqua
set output
exit ;
EOFMarker
fi

<< COMMENTOUT
# points/point[i].dat アニメーションの描画
set term gif animate optimize delay 4 size 480,360
set output '$dirname2/ms.gif'
set nokey
set xrange[-1:$N]
set yrange[-1:$N]
dt = 0.0025
t = 0
n = 0
skip = 50
do for [n = 0:(11000/skip-1)] {
  t = n*dt*skip
  set label 1 center at screen 0.5,0.9 sprintf("t=%.4f",t)
  # plot n th row in all i th file ( in the same frame )
  plot for [i=0:$L-1] '$dirname1/data/points/point'.i .'.dat' every ::((n+1)*skip)::((n+1)*skip) with points pt 7 lc 'dark-violet'
#  if( n%100 == 1 ) print("complete time step %d\n",n+1);
}
set out
#
set terminal aqua
set output
exit ;
EOFMarker
fi
COMMENTOUT
