# 行列Lの確認プロット
M=12
i=0
plot for[i=2:*] "L.dat" using ($1*0.0025):i w l title"l_{".(i-2)."}"
