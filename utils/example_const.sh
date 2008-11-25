./constDist X.dat 0.0 0.0 > constDistr_out.dat

echo "set term png; set output \"example_const.png\"; plot \"constDistr_out.dat\" w lines"|gnuplot
