#set terminal png
#set output "outputs.png"
plot for [i=2:4] "outputs.dat" using ($1*0.0025):i w l title sprintf('ms-%d',i-1)
replot for [i=5:7] "outputs.dat" using ($1*0.0025):i w l title sprintf('narma-%d',i-4)
#set terminal aqua
#set output
