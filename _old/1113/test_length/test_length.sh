cd `dirname $0` # change current directry to [ms]
gcc test_length.c -o test_length
echo comple completed.
mkdir -p ./points
mkdir -p ./springs
mkdir -p ./waves
echo simulating mass-spring system...
./test_length #./ms N NSTEP
echo plotting waves...
gnuplot -persist <<-EOFMarker
  load "test_length.plt" ;
  exit ;
EOFMarker
