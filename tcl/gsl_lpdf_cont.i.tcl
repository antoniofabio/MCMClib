##RUN THIS TCL SCRIPT TO GENERATE THE CORRESPONDING SWIG '.i' FILE

proc TYPE_PAR {prefix} {
	set tmp _lpdf_p
	return $prefix$tmp
}

proc DECLARE_2PAR {prefix par1 par2} {
	set tmp_lpdf _lpdf
	set tmp_alloc _lpdf_alloc
	set tmp_free _lpdf_free
	set ans {
typedef struct {
	double* $par1;
	double* $par2;
} [TYPE_PAR $prefix];

[TYPE_PAR $prefix]* mcmclib_$prefix$tmp_alloc (double* $par1, double* $par2);
void mcmclib_$prefix$tmp_free ([TYPE_PAR $prefix]* p);
%constant double mcmclib_$prefix$tmp_lpdf (gsl_vector* x, void* in_p);
}
	puts [subst $ans]
}

proc DECLARE_1PAR {prefix par1} {
	set tmp_lpdf _lpdf
	set tmp_alloc _lpdf_alloc
	set tmp_free _lpdf_free
	set ans {
typedef struct {
	double* $par1;
} [TYPE_PAR $prefix];

[TYPE_PAR $prefix]* mcmclib_$prefix$tmp_alloc (double* $par1);
void mcmclib_$prefix$tmp_free ([TYPE_PAR $prefix]* p);
%constant double mcmclib_$prefix$tmp_lpdf (gsl_vector* x, void* in_p);
}
	puts [subst $ans]
}

DECLARE_1PAR gaussian  sd
DECLARE_1PAR exponential  mean
DECLARE_1PAR laplace  a
DECLARE_2PAR exppow  a  b
DECLARE_1PAR cauchy  a
DECLARE_1PAR rayleigh  sigma
DECLARE_2PAR rayleigh_tail  a  sigma
DECLARE_2PAR gamma  a  b
DECLARE_2PAR flat  a  b
DECLARE_2PAR lognormal  zeta  sigma
DECLARE_1PAR chisq  nu
DECLARE_2PAR fdist  nu1  nu2
DECLARE_1PAR tdist  nu
DECLARE_2PAR beta  a  b
DECLARE_1PAR logistic  a
DECLARE_2PAR pareto  a  b
DECLARE_2PAR weibull  a  b
DECLARE_2PAR gumbel1  a  b
DECLARE_2PAR gumbel2  a  b
