load ./mcmclib.so
source ./util.tcl

##vector_list
set lst [mcmclib_vector_list_alloc]
stopifne [vector_list_str_next_get $lst] NULL
set v [l2v [list 1 2 3]]
set last [mcmclib_vector_list_append $v $lst]
stopifne [vector_list_str_v_get $last] $v

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
