cd `dirname $0` # change current directry to [ms]
gcc -I/usr/local/include -c ms.c -o ms.o
gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#echo comple completed.
# [0]./ms, [1]debug_flag [2]N(pow of integer), [3]mu_k, [4]sigma_k, [5]mu_g, [6]sigma_g, [7] dirname
# 正規分布の muは平均、sigmaは標準偏差
#time=`date "+%m%d_%H%M%S"`
#dirname1="./${time}"
dirname1="testSample"
N=6
debug=1
mu_k=100
sigma_k=30
mu_g=0.1
sigma_g=0.0
d_idx=0
L=$(($N**2)) #点の個数
M=$(($(($(($N-1))*$N))*2)) # ばねの本数
mkdir -p ${dirname1}

if [ $debug -eq 0 ] ; then
echo "# N=$N, M=$M" > ${dirname1}/data.txt
echo '# mu_k, sigma_k, mu_g, sigma_g, le_avr, err_narma2, err_narma10, err_narma20, err_narma30' >> ${dirname1}/data.txt
for i in `seq 1 100`
do
mu_g=`echo "scale=5; 5.0 / 100.0" | bc` #scale=小数点以下の桁数
mu_g=`echo "scale=5; $mu_g*$i" | bc`
sigma_g=`echo "scale=5; 1.0 / 100.0" | bc`
sigma_g=`echo "scale=5; $sigma_g*$i" | bc`
#mu_k=$((10*$i))
#sigma_k=$((3*$i))
for j in `seq 1 2` #同じパラメータで2回ずつ繰り返し
do
echo "mu_k=$mu_k, sigma_k=$sigma_k"
./ms $debug $N $mu_k $sigma_k $mu_g $sigma_g ${dirname1}
done
done

gnuplot -persist <<-EOFMarker
set terminal png
set xlabel 'lyapunov exponent'
set ylabel 'err'
set yrange [0:0.0001]
set output 'k/le-err-2.png'
plot 'k/data.txt' using 5:6 w p title 'NARMA2'
set yrange [0:0.1]
set output 'k/le-err-10.png'
plot 'k/data.txt' using 5:7 w p title 'NARMA10'
set yrange [0:0.1]
set output 'k/le-err-20.png'
plot 'k/data.txt' using 5:8 w p title 'NARMA20'
set yrange [0:10]
set output 'k/le-err-30.png'
plot 'k/data.txt' using 5:9 w p title 'NARMA30'
#
set terminal aqua
set output
exit ;
EOFMarker
fi

if [ $debug -eq 1 ] ; then
mkdir -p ${dirname1}/data #outputs, le .dat
mkdir -p ${dirname1}/data/points
mkdir -p ${dirname1}/data/springs
mkdir -p ${dirname1}/results # parameters and results
mkdir -p ${dirname1}/results/pictures
mkdir -p ${dirname1}/results/pictures/springs
#mkdir -p ${dirname1}/results/pictures/convergence
dirname2="${dirname1}/results/pictures"
./ms $debug $N $mu_k $sigma_k $mu_g $sigma_g ${dirname1}

gnuplot -persist <<-EOFMarker
set terminal png
# data/springs/length[i].dat のプロット
i=0
set xlabel 'real time'
set xrange [52.5:57.5]
do for [i=0:$M-1]{
set output sprintf('$dirname2/springs/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('$dirname1/data/springs/length%d.dat',i) using 2:3 w lp title sprintf('l_{%d}(t)',i)
}
set ylabel 'bias'
set output '$dirname2/springs/bais.png'
plot '$dirname1/data/springs/bias.dat' using 2:3 w l notitle

# data/input.dat 入力時系列のプロット
set ylabel 'I(t)'
set xrange [52.5:57.5] #1000step
#set xrange [6000:8000]
set output '$dirname2/input.png'
plot '$dirname1/data/input.dat' using 2:3 w lp notitle
unset xrange

# data/outputs.dat NARMAモデルと近似波形のプロット
set xrange [50:105]
set ylabel 'y(t)'
set output '$dirname2/NARMA2.png'
plot '$dirname1/data/outputs.dat' using 2:3 w lp title "ms-NARMA2", '$dirname1/data/outputs.dat' using 2:7 w lp title "NARMA2"
set output '$dirname2/NARMA10.png'
plot '$dirname1/data/outputs.dat' using 2:4 w lp title "ms-NARMA10", '$dirname1/data/outputs.dat' using 2:8 w lp title "NARMA10"
set output '$dirname2/NARMA20.png'
plot '$dirname1/data/outputs.dat' using 2:5 w lp title "ms-NARMA20", '$dirname1/data/outputs.dat' using 2:9 w lp title "NARMA20"
set output '$dirname2/NARMA30.png'
plot '$dirname1/data/outputs.dat' using 2:6 w lp title "ms-NARMA30", '$dirname1/data/outputs.dat' using 2:10 w lp title "NARMA30"

set xrange [52.5:57.5] #1000step
set output '$dirname2/NARMA2_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:3 w lp title "ms-NARMA2", '$dirname1/data/outputs.dat' using 2:7 w lp title "NARMA2"
set output '$dirname2/NARMA10_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:4 w lp title "ms-NARMA10", '$dirname1/data/outputs.dat' using 2:8 w lp title "NARMA10"
set output '$dirname2/NARMA20_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:5 w lp title "ms-NARMA20", '$dirname1/data/outputs.dat' using 2:9 w lp title "NARMA20"
set output '$dirname2/NARMA30_zoom.png'
plot '$dirname1/data/outputs.dat' using 2:6 w lp title "ms-NARMA30", '$dirname1/data/outputs.dat' using 2:10 w lp title "NARMA30"
unset xrange

# results/le.dat　のプロット
set ylabel 'lyapunov exponent'
set output '$dirname2/le.png'
plot '$dirname1/data/le.dat' using 1:3 w l notitle

#
set terminal aqua
set output
exit ;
EOFMarker
fi

<< COMMENTOUT
# data/springs/length[i].dat のプロット(ずらさない系lとずらした系l_dの比較用)
set xrange [0:2.5]
do for [i=0:$M-1]{
set out sprintf('$dirname2/convergence/l_%03d.png',i)
set ylabel sprintf('l_{%d}(t)',i)
plot sprintf('$dirname1/data/springs/length%d.dat',i) using 2:3 w l title sprintf('l_{%d}(t)',i), sprintf('$dirname1/data/springs/length%d.dat',i) using 2:4 w l title sprintf('dif-l_{%d}(t) (index:$d_idx)',i)
}

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
COMMENTOUT
