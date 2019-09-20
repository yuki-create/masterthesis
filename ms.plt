set term gif animate optimize delay 4 size 480,360
set output './data25/ms.gif'
set nokey
set xrange[-10:20]
set yrange[-10:20]
dt = 0.01
t = 0
n = 0
NSTEP = 1000
N = 25

do for [n = 0:NSTEP ] {
  t = n*dt
  set label 1 center at screen 0.5,0.9 sprintf("t=%.2f",t)
  # plot n th row in all i th file ( in the same frame )
  plot for [i=0:N-1] './data25/data_point'.i .'.dat' every ::(n+1)::(n+1)
#  if( n%100 == 1 ) print("complete time step %d\n",n+1);
}
set out
