dirname1="./test"
N=3
debug=1
L=$(($N**2)) #点の個数
M=$(($(($(($N-1))*$N))*2)) # ばねの本数
end=$((L*L-1))
dirname2="${dirname1}/pictures"
tmp=20
#gcc -I/usr/local/include -c ms.c -o ms.o
#gcc -L/Users/saki/lapack-3.8.0 ms.o -llapacke -llapack -lcblas -lblas -lgfortran -lm -o ms
#./ms $debug $N 800 100 0.03 0.01 ${dirname1}

IFS=$'\n' #区切り文字をここで変更している
list=(`cat ${dirname1}/results/graph.dat`)
#i=0
#for i in ${list[@]}
#do
 #echo ${list[i]}
#done

echo "COUNT=${#list[@]}"
gnuplot -persist <<-EOFMarker
i=0
j=0
# データサイズの取得
#stats "$dirname1/results/graph.dat" nooutput
#r = STATS_records # データの行数
#c = STATS_columns # データの列数（コラム数）
# 配列の宣言
array G[$L*$L]
# データの配列への保存
#stats "$dirname1/results/graph.dat" using (sum[i=1:c] (A_2d[i + $0*c] = column(i), 0)) nooutput
# 配列の内容の表示
#do for [i=0:$L-1]{
#  do for [j=0:$L-1]{
#    idx=j+i*$L
#      G[j+i*$L+1] = ${list[j+i*$L]}
#     print idx
#     print ${list[idx]}
#     }
#  }
#do for [i=0:$end]{
#  print i,`echo ${list[i]}`,${list[1]}
#}
#print `sed -n 10p "${dirname1}/results/graph.dat"`
i=0
print $tmp
`$tmp=`
print $tmp
#do for [i=0:$end]{
#  print "$dirname1/results/graph.dat"  every ::i::i using 1
#}
exit ;
EOFMarker
<< COMMENTOUT
gnuplot -persist <<-EOFMarker
i=0
# points/point[i].dat アニメーションの描画
set term gif animate optimize delay 4 size 480,360
set output '$dirname2/ms.gif'
set nokey
set xrange[-1:$N]
set yrange[-1:$N]
dt = 0.0025
t = 0
n = 0
skip = 50
do for [n = 0:(11000/skip-1)] {
  t = n*dt*skip
  set label 1 center at screen 0.5,0.9 sprintf("t=%.4f",t)
  # plot n th row in all i th file ( in the same frame )
  plot for [i=0:$L-1] '$dirname1/points/point'.i .'.dat' every ::((n+1)*skip)::((n+1)*skip)
#  if( n%100 == 1 ) print("complete time step %d\n",n+1);
}
set out
#
set terminal aqua
set output
exit ;
EOFMarker
COMMENTOUT
