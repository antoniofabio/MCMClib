load ./mcmclib.so
source ./util.tcl

set pg [mcmclib_gaussian_lpdf_alloc [ll2a 1.0]]
set x [l2v 0.0]
mcmclib_gaussian_lpdf $x $pg
set cb $mcmclib_gaussian_lpdf_cb
mcmclib_gaussian_lpdf_free $pg
