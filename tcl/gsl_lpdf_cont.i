
typedef struct {
	double* sd;
} gaussian_lpdf_p;

gaussian_lpdf_p* mcmclib_gaussian_lpdf_alloc (double* sd);
void mcmclib_gaussian_lpdf_free (gaussian_lpdf_p* p);
%callback("%s_cb");
double mcmclib_gaussian_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* mean;
} exponential_lpdf_p;

exponential_lpdf_p* mcmclib_exponential_lpdf_alloc (double* mean);
void mcmclib_exponential_lpdf_free (exponential_lpdf_p* p);
%callback("%s_cb");
double mcmclib_exponential_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
} laplace_lpdf_p;

laplace_lpdf_p* mcmclib_laplace_lpdf_alloc (double* a);
void mcmclib_laplace_lpdf_free (laplace_lpdf_p* p);
%callback("%s_cb");
double mcmclib_laplace_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} exppow_lpdf_p;

exppow_lpdf_p* mcmclib_exppow_lpdf_alloc (double* a, double* b);
void mcmclib_exppow_lpdf_free (exppow_lpdf_p* p);
%callback("%s_cb");
double mcmclib_exppow_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
} cauchy_lpdf_p;

cauchy_lpdf_p* mcmclib_cauchy_lpdf_alloc (double* a);
void mcmclib_cauchy_lpdf_free (cauchy_lpdf_p* p);
%callback("%s_cb");
double mcmclib_cauchy_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* sigma;
} rayleigh_lpdf_p;

rayleigh_lpdf_p* mcmclib_rayleigh_lpdf_alloc (double* sigma);
void mcmclib_rayleigh_lpdf_free (rayleigh_lpdf_p* p);
%callback("%s_cb");
double mcmclib_rayleigh_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* sigma;
} rayleigh_tail_lpdf_p;

rayleigh_tail_lpdf_p* mcmclib_rayleigh_tail_lpdf_alloc (double* a, double* sigma);
void mcmclib_rayleigh_tail_lpdf_free (rayleigh_tail_lpdf_p* p);
%callback("%s_cb");
double mcmclib_rayleigh_tail_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} gamma_lpdf_p;

gamma_lpdf_p* mcmclib_gamma_lpdf_alloc (double* a, double* b);
void mcmclib_gamma_lpdf_free (gamma_lpdf_p* p);
%callback("%s_cb");
double mcmclib_gamma_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} flat_lpdf_p;

flat_lpdf_p* mcmclib_flat_lpdf_alloc (double* a, double* b);
void mcmclib_flat_lpdf_free (flat_lpdf_p* p);
%callback("%s_cb");
double mcmclib_flat_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* zeta;
	double* sigma;
} lognormal_lpdf_p;

lognormal_lpdf_p* mcmclib_lognormal_lpdf_alloc (double* zeta, double* sigma);
void mcmclib_lognormal_lpdf_free (lognormal_lpdf_p* p);
%callback("%s_cb");
double mcmclib_lognormal_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* nu;
} chisq_lpdf_p;

chisq_lpdf_p* mcmclib_chisq_lpdf_alloc (double* nu);
void mcmclib_chisq_lpdf_free (chisq_lpdf_p* p);
%callback("%s_cb");
double mcmclib_chisq_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* nu1;
	double* nu2;
} fdist_lpdf_p;

fdist_lpdf_p* mcmclib_fdist_lpdf_alloc (double* nu1, double* nu2);
void mcmclib_fdist_lpdf_free (fdist_lpdf_p* p);
%callback("%s_cb");
double mcmclib_fdist_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* nu;
} tdist_lpdf_p;

tdist_lpdf_p* mcmclib_tdist_lpdf_alloc (double* nu);
void mcmclib_tdist_lpdf_free (tdist_lpdf_p* p);
%callback("%s_cb");
double mcmclib_tdist_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} beta_lpdf_p;

beta_lpdf_p* mcmclib_beta_lpdf_alloc (double* a, double* b);
void mcmclib_beta_lpdf_free (beta_lpdf_p* p);
%callback("%s_cb");
double mcmclib_beta_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
} logistic_lpdf_p;

logistic_lpdf_p* mcmclib_logistic_lpdf_alloc (double* a);
void mcmclib_logistic_lpdf_free (logistic_lpdf_p* p);
%callback("%s_cb");
double mcmclib_logistic_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} pareto_lpdf_p;

pareto_lpdf_p* mcmclib_pareto_lpdf_alloc (double* a, double* b);
void mcmclib_pareto_lpdf_free (pareto_lpdf_p* p);
%callback("%s_cb");
double mcmclib_pareto_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} weibull_lpdf_p;

weibull_lpdf_p* mcmclib_weibull_lpdf_alloc (double* a, double* b);
void mcmclib_weibull_lpdf_free (weibull_lpdf_p* p);
%callback("%s_cb");
double mcmclib_weibull_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} gumbel1_lpdf_p;

gumbel1_lpdf_p* mcmclib_gumbel1_lpdf_alloc (double* a, double* b);
void mcmclib_gumbel1_lpdf_free (gumbel1_lpdf_p* p);
%callback("%s_cb");
double mcmclib_gumbel1_lpdf (gsl_vector* x, void* in_p);
%nocallback;


typedef struct {
	double* a;
	double* b;
} gumbel2_lpdf_p;

gumbel2_lpdf_p* mcmclib_gumbel2_lpdf_alloc (double* a, double* b);
void mcmclib_gumbel2_lpdf_free (gumbel2_lpdf_p* p);
%callback("%s_cb");
double mcmclib_gumbel2_lpdf (gsl_vector* x, void* in_p);
%nocallback;

