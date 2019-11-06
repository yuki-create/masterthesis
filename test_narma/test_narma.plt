set terminal png
set output "input.png"
plot "outputs.dat" using($1*0.0025):2 title"input,dt=0.0025,T=1,NSTEP=10000" with lines
set output "narma2.png"
set yrange [0.18:0.2]
plot "outputs.dat" using($1*0.0025):3 title"narma2,dt=0.0025,T=1,NSTEP=10000" with lines
set output "narma10.png"
set yrange [0.14:0.26]
plot "outputs.dat" using($1*0.0025):4 title"narma10,dt=0.0025,T=1,NSTEP=10000" with lines
set output "narma20.png"
set yrange [0.16:0.26]
plot "outputs.dat" using($1*0.0025):5 title"narma20,dt=0.0025,T=1,NSTEP=10000" with lines

unset yrange
set xrange [12.5:13.75]
set output "input_t[12.5:13.75].png"
plot "outputs.dat" using($1*0.0025):2 title"input,dt=0.0025,T=1,NSTEP=10000" with lines
set output "narma2_t[12.5:13.75].png"
set yrange [0.18:0.2]
plot "outputs.dat" using($1*0.0025):3 title"narma2,dt=0.0025,T=1,NSTEP=10000" with lines
set output "narma10_t[12.5:13.75].png"
set yrange [0.14:0.26]
plot "outputs.dat" using($1*0.0025):4 title"narma10,dt=0.0025,T=1,NSTEP=10000" with lines
set output "narma20_t[12.5:13.75].png"
set yrange [0.17:0.24]
plot "outputs.dat" using($1*0.0025):5 title"narma20,dt=0.0025,T=1,NSTEP=10000" with lines

unset yrange
unset xrange
set terminal aqua
set output
