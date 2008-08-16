load ./mcmclib.so
source ./util.tcl

set pg [mcmclib_gaussian_lpdf_alloc [l2a 1.0]]
set x [l2v 0.0]
set value [mcmclib_gaussian_lpdf $x $pg]
stopifne [expr round($value * 1e6)] -918939
set cb $mcmclib_gaussian_lpdf_cb
mcmclib_gaussian_lpdf_free $pg
