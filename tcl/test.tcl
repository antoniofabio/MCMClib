load ./mcmclib.so
source ./util.tcl

##vector_list
set ll [list\
	[list 1.0 2.0 3.0]\
	[list 2.0 4.0 6.0]\
	[list 3.0 6.0 9.0]\
	[list 4.0 8.0 12.0]]
set vl [ll2vl $ll]
stopifne [mcmclib_vector_list_length $vl] 4
stopifne [vl2ll $vl] $ll

##list<->matrix
set ll {{1.0 2.0 3.0} {2.0 4.0 6.0} {3.0 6.0 9.0} {4.0 8.0 12.0}}
set m [ll2m $ll]
stopifne [m2ll $m] $ll

##init RNG
set r [gsl_rng_alloc $gsl_rng_default]

##gaussian random walk
set x [l2v 1]
mcmclib_gauss_rw $r $mcmclib_test_dunif $x NULL 1.0

set x [l2v 0.7]
set chain [list]
for {set i 0} {$i < 10000} {incr i} {
	mcmclib_gauss_rw $r $mcmclib_test_dunif $x NULL 0.1
	set lx [v2l $x]
	stopifnot [expr ($lx >= 0) & ($lx <= 1)]
	lappend chain $lx
}
