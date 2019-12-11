plot for [i=2:4] "test_outputsMS.dat" using ($1*0.0025):i w l title sprintf('ms-%d',i-1)
replot for [i=5:7] "test_outputsMS.dat" using ($1*0.0025):i w l title sprintf('narma-%d',i-4)
