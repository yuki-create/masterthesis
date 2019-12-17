cd `dirname $0` # change current directry to [ms]
#gcc ms.c -o test_length
gcc -I/usr/local/include -c exe.c -o exe.o
gcc -L/Users/saki/lapack-3.8.0 exe.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o exe
#echo comple completed.
str1="."
mkdir -p ${str1}/points
mkdir -p ${str1}/springs
mkdir -p ${str1}/results # outputs.dat, parameters and results
for i in `seq 1 1000`
do
  echo "$i 回目のループです。"
  ./exe
done
echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "./le_err.plt" ;
  exit ;
EOFMarker
