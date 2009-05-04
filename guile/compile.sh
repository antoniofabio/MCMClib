swig -guile mcmclib.i
gcc -shared --std=gnu99 ${CFLAGS} -fPIC mcmclib_wrap.c -I../src -lgsl -lgslcblas -lm ${LDFLAGS} ../src/libmcmclib.a -o libschememcmclib.so
