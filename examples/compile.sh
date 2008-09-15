gcc --std=gnu99 example1.c ../src/libmcmclib.a -I../src -lgsl -lgslcblas -o example1
gcc --std=gnu99 -g example2.c ../src/libmcmclib.a -I../src -lgsl -lgslcblas -o example2
