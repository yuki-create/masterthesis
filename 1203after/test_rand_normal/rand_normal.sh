gcc -I/usr/local/include -c rand_normal.c -o rand_normal.o
gcc -L/Users/saki/lapack-3.8.0 rand_normal.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o rand_normal
./rand_normal 9 0 1 0 1 | sort -n | uniq -c > rand_normal.txt
gnuplot -persist <<-EOFMarker
  plot "rand_normal.txt" using 2:1 with histeps
  exit ;
EOFMarker
