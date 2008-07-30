load ./mcmclib.so
source ./util.tcl

##vector_list
set lst [mcmclib_vector_list_alloc]
stopifne [vector_list_str_next_get $lst] NULL
set v [list2vector [list 1 2 3]]
mcmclib_vector_list_add $v $lst
stopifne [vector_list_str_v_get [vector_list_str_next_get $lst]] $v

##init RNG
set r [gsl_rng_alloc $gsl_rng_default]

##gaussian random walk
set x [list2vector 1]
mcmclib_gauss_rw $r $mcmclib_test_dunif $x NULL 1.0

set x [list2vector 0.7]
set chain [list]
for {set i 0} {$i < 10000} {incr i} {
	mcmclib_gauss_rw $r $mcmclib_test_dunif $x NULL 0.1
	set lx [vector2list $x]
	stopifnot [expr ($lx >= 0) & ($lx <= 1)]
	lappend chain $lx
}
