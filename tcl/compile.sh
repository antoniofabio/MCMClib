swig -w454 -tcl mcmclib.i
gcc -shared mcmclib_wrap.c -I../src -I/usr/include/tcl8.4 `pkg-config gsl --libs` ../src/libmcmclib.a -o mcmclib.so

