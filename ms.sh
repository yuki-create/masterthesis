cd `dirname $0` # change current directry to [ms]
gcc ms.c -o test_length
str1="test_length"
mkdir -p ./${str1}/points
mkdir -p ./${str1}/springs
./test #./ms N NSTEP
<< COMMENTOUT
echo drowing animation...
gnuplot -persist <<-EOFMarker
  load "ms.plt" ;
  exit ;
EOFMarker
COMMENTOUT
