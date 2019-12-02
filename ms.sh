cd `dirname $0` # change current directry to [ms]
#gcc ms.c -o test_length
gcc -I/usr/local/include -c ms.c -o ms.o
gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#echo comple completed.
# [0]./ms, [1]repeat [2]N(pow of integer), [3]mu_k, [4]sigma_k, [5]mu_g, [6]sigma_g
# 正規分布の muは平均、sigmaは標準偏差
str1="."
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
mkdir -p ${str1}/results # outputs.dat, parameters and results
echo -n > ${str1}/results/le.dat
./ms 10 9 800 100 0.01 0.01
gnuplot -persist <<-EOFMarker
  plot "./results/le.dat" w l ;
  exit ;
EOFMarker
<< COMMENTOUT
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
