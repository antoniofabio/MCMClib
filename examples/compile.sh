gcc --std=gnu99 -O0 -g $CFLAGS example1.c ../src/libmcmclib.a -I../src $LDFLAGS -lgsl -lgslcblas -lm -o example1
gcc --std=gnu99 -O0 -g $CFLAGS example2.c ../src/libmcmclib.a -I../src $LDFLAGS -lgsl -lgslcblas -lm -o example2
gcc --std=gnu99 -O0 -g $CFLAGS example3.c ../src/libmcmclib.a -I../src $LDFLAGS -lgsl -lgslcblas -lm -o example3

gcc --std=gnu99 -O0 -g $CFLAGS exRapt.c ../src/libmcmclib.a -I../src $LDFLAGS -lgsl -lgslcblas -lm -o exRapt
gcc --std=gnu99 -O0 -g $CFLAGS exRapt2.c ../src/libmcmclib.a -I../src $LDFLAGS -lgsl -lgslcblas -lm -o exRapt2

