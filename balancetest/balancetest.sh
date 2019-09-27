gcc balancetest.c -o balancetest
./balancetest #./ms N NSTEP
echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "balancetest.plt" ;
  exit ;
EOFMarker
