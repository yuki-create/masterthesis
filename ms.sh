cd `dirname $0` # change current directry to [ms]
#gcc ms.c -o test_length
gcc -I/usr/local/include -c ms.c -o ms.o
gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#echo comple completed.
# [0]./ms, [1]debug_frag [2]N(pow of integer), [3]mu_k, [4]sigma_k, [5]mu_g, [6]sigma_g, [7] dirname
# 正規分布の muは平均、sigmaは標準偏差
time=`date "+%m%d_%H%M%S"`
# str1="./${time}"
str1="./test"
mkdir ${str1}
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
mkdir -p ${str1}/results # outputs.dat, parameters and results
./ms 1 3 800 100 0.03 0.01 ${str1}
<< COMMENTOUT
gnuplot -persist <<-EOFMarker
  set terminal x11 1
  plot "./log1.dat" using 1:2 w l ;
  replot "./log1.dat" using 1:5 w l ;
  set terminal x11 2
  plot "./log1.dat" using 1:3 w l ;
  replot "./log1.dat" using 1:6 w l ;
  set terminal x11 3
  plot "./log1.dat" using 1:4 w l ;
  replot "./log1.dat" using 1:7 w l ;
  exit ;
EOFMarker

gnuplot -persist <<-EOFMarker
  plot "./log1.dat" using 1:2 w l ;
  replot "./log1.dat" using 1:3 w l ;
  exit ;
EOFMarker
gnuplot -persist <<-EOFMarker
  set xrange [100:11000] ;
  plot "./results/le.dat" w l ;
  exit ;
EOFMarker


gnuplot -persist <<-EOFMarker
  plot "test_length.txt" using 1:2 w l ;
  replot "test_length.txt" using 1:3 w l ;
  exit ;
EOFMarker

gnuplot -persist <<-EOFMarker
  set xrange [4000:5500]
  plot "test_force.txt" using 1:2 w l;
  replot "test_force.txt" using 1:3 w l;
  exit ;
EOFMarker

str1="."
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
mkdir -p ${str1}/results # outputs.dat, parameters and results
#echo simulating ms...
./ms #./ms N NSTEP

echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "./outputs.plt" ;
  exit ;
EOFMarker

# plot lyapunov exponent
echo drowing animation...
gnuplot -persist <<-EOFMarker
  set xrange [100:11000] ;
  plot "./results/le.dat" with l
  exit ;
EOFMarker
COMMENTOUT
