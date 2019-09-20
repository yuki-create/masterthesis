gcc ms.c -o ms
./ms #./ms N NSTEP
echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "ms.plt" ;
  exit ;
EOFMarker
