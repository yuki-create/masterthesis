cd `dirname $0` # change current directry to [test_outputsMS]
gcc -I/usr/local/include -c test_outputsMS.c -o test_outputsMS.o
gcc -L/Users/saki/lapack-3.8.0 test_outputsMS.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o test_outputsMS
#echo comple completed.
str1="."
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
#echo simulating test_outputsMS...
./test_outputsMS > test_outputsMS.dat
#echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "test_outputsMS.plt" ;
  exit ;
EOFMarker
