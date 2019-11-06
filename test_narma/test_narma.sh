# msディレクトリ下で ./test_narma/test_narma.sh を実行
gcc ms.c -o ./test_narma/test_narma
./test_narma/test_narma > ./test_narma/test_narma.dat
gnuplot -persist <<-EOFMarker
  load "./test_narma/test_narma.plt" ;
  exit ;
EOFMarker
