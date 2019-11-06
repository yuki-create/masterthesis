cd `dirname $0` # change current directry to [test_narma]
gcc test_narma.c -o test_narma
mkdir -p ./points
mkdir -p ./springs
./test_narma > ./outputs.dat
gnuplot -persist <<-EOFMarker
  load "test_narma.plt" ;
  exit ;
EOFMarker
