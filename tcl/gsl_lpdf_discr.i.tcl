##RUN THIS TCL SCRIPT TO GENERATE THE CORRESPONDING SWIG '.i' FILE

proc TYPE_PAR {prefix} {
	set tmp _lpdf_p
	return $prefix$tmp
}

proc DECLARE_2PAR {prefix type1 par1 type2 par2} {
	set tmp_lpdf _lpdf
	set tmp_alloc _lpdf_alloc
	set tmp_free _lpdf_free
	set ans {
typedef struct {
	$type1 * $par1;
	$type2 * $par2;
} [TYPE_PAR $prefix];

[TYPE_PAR $prefix]* mcmclib_$prefix$tmp_alloc ($type1 * $par1, $type2 * $par2);
void mcmclib_$prefix$tmp_free ([TYPE_PAR $prefix]* p);
%callback("%s_cb");
double mcmclib_$prefix$tmp_lpdf (gsl_vector* x, void* in_p);
%nocallback;
}
	puts [subst $ans]
}

proc DECLARE_1PAR {prefix type1 par1} {
	set tmp_lpdf _lpdf
	set tmp_alloc _lpdf_alloc
	set tmp_free _lpdf_free
	set ans {
typedef struct {
	$type1 * $par1;
} [TYPE_PAR $prefix];

[TYPE_PAR $prefix]* mcmclib_$prefix$tmp_alloc ($type1 * $par1);
void mcmclib_$prefix$tmp_free ([TYPE_PAR $prefix]* p);
%callback("%s_cb");
double mcmclib_$prefix$tmp_lpdf (gsl_vector* x, void* in_p);
%nocallback;
}
	puts [subst $ans]
}

DECLARE_1PAR poisson double mu
DECLARE_1PAR bernoulli double p
DECLARE_2PAR binomial double p int n
DECLARE_2PAR negative_binomial double p int n
DECLARE_2PAR pascal double p int n
DECLARE_1PAR geometric double p
DECLARE_1PAR logarithmic double p
