./distrDist 2 X.dat 1000 Y.dat 10000 X0.dat 50 > distrDist_out.dat

echo "set term png; set output \"distrDist_example.png\"; plot \"distrDist_out.dat\" w lines"|gnuplot
