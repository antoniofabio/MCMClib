gcc -O0 -g --std=gnu99 $CFLAGS constDist.c $LDFLAGS -lgsl -lgslcblas -lm -o constDist
gcc -O0 -g --std=gnu99 $CFLAGS distrDist.c $LDFLAGS -lgsl -lgslcblas -lm -o distrDist

