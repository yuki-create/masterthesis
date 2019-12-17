plot for [i=2:4] "./results/outputs.dat" using ($1*0.0025):i w l title sprintf('ms-%d',i-1)
replot for [i=5:7] "./results/outputs.dat" using ($1*0.0025):i w l title sprintf('narma-%d',i-4)
set terminal png
set out "./results/outputs.png"
replot
set terminal aqua
set output
