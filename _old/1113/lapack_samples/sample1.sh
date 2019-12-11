gcc -I/usr/local/include -c sample1.c -o sample1.o
gcc -L/Users/saki/lapack-3.8.0 sample1.o  -llapacke -llapack -lcblas -lblas -lgfortran -lm -o sample1
./sample1
