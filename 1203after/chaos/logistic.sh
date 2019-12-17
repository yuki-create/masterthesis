gcc logistic.c -o logistic
./logistic > logistic_le.dat
gnuplot -persist <<-EOFMarker
  set xlabel "n"
  set ylabel "le"
  plot "./logistic_le.dat" using 1:2 with l title"a=4.0,lyapunov"
  replot "./logistic_le.dat" using 1:3 with l title"a=4.0,accurate le"
  exit ;
EOFMarker
