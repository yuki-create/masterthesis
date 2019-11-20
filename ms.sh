cd `dirname $0` # change current directry to [ms]
#gcc ms.c -o test_length
gcc -I/usr/local/include -c ms.c -o ms.o
gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#echo comple completed.
str1="."
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
mkdir -p ${str1}/results # outputs.dat, parameters and results
#echo simulating ms...
./ms #./ms N NSTEP
<< COMMENTOUT
echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "./outputs.plt" ;
  exit ;
EOFMarker
COMMENTOUT
# plot lyapunov exponent
echo drowing animation...
gnuplot -persist <<-EOFMarker
  set xrange [1000:16000] ;
  plot "./results/le.dat" with l
  exit ;
EOFMarker
