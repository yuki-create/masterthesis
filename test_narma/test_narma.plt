set terminal png
set output "./test_narma/input.png"
plot "./test_narma/test_narma.dat" using($1*0.01):2 title"input,dt=0.0025,T=1,NSTEP=10000" with lines
set output "./test_narma/narma2.png"
set yrange [0.18:0.2]
plot "./test_narma/test_narma.dat" using($1*0.01):3 title"narma2,dt=0.0025,T=1,NSTEP=10000" with lines
set output "./test_narma/narma10.png"
set yrange [0.14:0.26]
plot "./test_narma/test_narma.dat" using($1*0.01):4 title"narma10,dt=0.0025,T=1,NSTEP=10000" with lines
set output "./test_narma/narma20.png"
set yrange [0.16:0.26]
plot "./test_narma/test_narma.dat" using($1*0.01):5 title"narma20,dt=0.0025,T=1,NSTEP=10000" with lines

unset yrange
set xrange [50:55]
set output "./test_narma/input_t50-55.png"
plot "./test_narma/test_narma.dat" using($1*0.01):2 title"input,dt=0.0025,T=1,NSTEP=10000" with lines
set output "./test_narma/narma2_t50-55.png"
set yrange [0.18:0.2]
plot "./test_narma/test_narma.dat" using($1*0.01):3 title"narma2,dt=0.0025,T=1,NSTEP=10000" with lines
set output "./test_narma/narma10_t50-55.png"
set yrange [0.14:0.26]
plot "./test_narma/test_narma.dat" using($1*0.01):4 title"narma10,dt=0.0025,T=1,NSTEP=10000" with lines
set output "./test_narma/narma20_t50-55.png"
set yrange [0.17:0.24]
plot "./test_narma/test_narma.dat" using($1*0.01):5 title"narma20,dt=0.0025,T=1,NSTEP=10000" with lines

unset yrange
unset xrange
set terminal aqua
set output
