gcc -I/usr/local/include -c matvec.c -o matvec.o
gcc -L/Users/saki/lapack-3.8.0 matvec.o  -llapacke -llapack -lcblas -lblas -lgfortran -lm -o matvec
./matvec
