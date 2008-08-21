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
stopifne [llength [vl2ll [mcmclib_vector_list_transpose $vl]]] 3
vl2ll [mcmclib_vector_list_transpose $vl]
stopifne [vl2ll [mcmclib_vector_list_transpose [mcmclib_vector_list_transpose $vl]]] [vl2ll $vl]

##list<->matrix
set ll {{1.0 2.0 3.0} {2.0 4.0 6.0} {3.0 6.0 9.0} {4.0 8.0 12.0}}
set m [ll2m $ll]
stopifne [m2ll $m] $ll

##variance/covariance matrix
set m [ll2m {{1.0 2.0 3.0} {2.0 4.0 6.0} {3.0 6.0 9.0} {4.0 8.0 12.0}}]
set mcov [gsl_matrix_alloc 3 3]
mcmclib_matrix_covariance $m $mcov
stopifne [m2ll $mcov] {{1.25 2.5 3.75} {2.5 5.0 7.5} {3.75 7.5 11.25}}

set m [ll2m {{1.0 2.0 3.0} {2.0 4.0 6.0} {3.0 6.0 9.0}}]
mcmclib_matrix_covariance $m $mcov
set n [new_intArray 1]
intArray_setitem $n 0 3
set mean [l2v "0 0 0"]
mcmclib_matrix_colmeans $m $mean
mcmclib_covariance_update $mcov $mean $n [l2v {4.0 8.0 12.0}]
stopifne [intArray_getitem $n 0] 4
stopifne [m2ll $mcov] {{1.25 2.5 3.75} {2.5 5.0 7.5} {3.75 7.5 11.25}}

##init RNG
set r [gsl_rng_alloc $gsl_rng_default]

##Multivariate Gaussian variates
set sigma [ll2m { {1 0} {0 1} }]
set out [gsl_vector_alloc 2]
set ans [list]
for {set i 0} {$i < 1000} {incr i} {
	mcmclib_mvnorm $r $sigma $out
	lappend ans [v2l $out]
}
set cov [gsl_matrix_alloc 2 2]
mcmclib_matrix_covariance [ll2m $ans] $cov
stopifne [m2ll $cov] \
	{{1.0020498344099165 0.059333707261511485} {0.059333707261511485 0.891687139422672}}

##gaussian random walk
set r [gsl_rng_alloc $gsl_rng_default]
set x [l2v 1]
set extra [mcmclib_gauss_rw_alloc 0.1 1]
mcmclib_gauss_rw $r $mcmclib_test_dunif $x NULL $extra

set x [l2v 0.7]
set chain [list]
for {set i 0} {$i < 10000} {incr i} {
	mcmclib_gauss_rw $r $mcmclib_test_dunif $x NULL $extra
	set lx [v2l $x]
	stopifnot [expr ($lx >= 0) & ($lx <= 1)]
	lappend chain $lx
}
mcmclib_gauss_rw_free $extra

##Adaptive Metropolis (Haario et al., 2001)
set r [gsl_rng_alloc $gsl_rng_default]
set x [l2v "0.5 0.5 0.5"]
set extra [mcmclib_gauss_am_alloc [ll2m { {0.1 0 0} {0 0.1 0} {0 0 0.1} }] 100]
mcmclib_gauss_am $r $mcmclib_test_dunif $x NULL $extra
set chain [list]
for {set i 0} {$i < 10000} {incr i} {
	mcmclib_gauss_am $r $mcmclib_test_dunif $x NULL $extra
	set lx [v2l $x]
	lappend chain $lx
}
set cov [gsl_matrix_alloc 3 3]
mcmclib_matrix_covariance [ll2m $chain] $cov
stopifne [m2ll $cov] \
	{{0.08415137716635926 -0.0014242134920127044 -5.698656336101036e-5}\
 {-0.0014242134920127044 0.08581935002944255 -0.0021620104951330468}\
 {-5.698656336101036e-5 -0.0021620104951330468 0.08401825041995266}}
