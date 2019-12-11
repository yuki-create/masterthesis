set terminal png
set xlabel "lyapunov exponent"
set ylabel "err"

set yrange [0:0.001]
set output './results/le-err-2.png'
plot "./results/le_err.dat" using 1:2 title "narma2"

set yrange [0:0.1]
set output './results/le-err-10.png'
plot "./results/le_err.dat" using 1:3 title "narma10"

set output './results/le-err-20.png'
plot "./results/le_err.dat" using 1:4 title "narma20"

set yrange [0:0.0001]
set output './results/zoom_le-err-2.png'
plot "./results/le_err.dat" using 1:2 title "narma2"

set yrange [0:0.03]
set output './results/zoom_le-err-10.png'
plot "./results/le_err.dat" using 1:3 title "narma10"

set yrange [0:0.02]
set output './results/zoom_le-err-20.png'
plot "./results/le_err.dat" using 1:4 title "narma20"

set terminal aqua
set output
