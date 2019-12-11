gcc -I/usr/local/include -c least_squares.c -o least_squares.o
gcc -L/Users/saki/lapack-3.8.0 least_squares.o  -llapacke -llapack -lcblas -lblas -lgfortran -lm -o least_squares
./least_squares
