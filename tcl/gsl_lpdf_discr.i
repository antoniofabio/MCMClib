%callback("%s_cb");

typedef struct {
	double * mu;
} poisson_lpdf_p;

poisson_lpdf_p* mcmclib_poisson_lpdf_alloc (double * mu);
void mcmclib_poisson_lpdf_free (poisson_lpdf_p* p);
%constant double mcmclib_poisson_lpdf (gsl_vector* x, void* in_p);


typedef struct {
	double * p;
} bernoulli_lpdf_p;

bernoulli_lpdf_p* mcmclib_bernoulli_lpdf_alloc (double * p);
void mcmclib_bernoulli_lpdf_free (bernoulli_lpdf_p* p);
%constant double mcmclib_bernoulli_lpdf (gsl_vector* x, void* in_p);


typedef struct {
	double * p;
	int * n;
} binomial_lpdf_p;

binomial_lpdf_p* mcmclib_binomial_lpdf_alloc (double * p, int * n);
void mcmclib_binomial_lpdf_free (binomial_lpdf_p* p);
%constant double mcmclib_binomial_lpdf (gsl_vector* x, void* in_p);


typedef struct {
	double * p;
	int * n;
} negative_binomial_lpdf_p;

negative_binomial_lpdf_p* mcmclib_negative_binomial_lpdf_alloc (double * p, int * n);
void mcmclib_negative_binomial_lpdf_free (negative_binomial_lpdf_p* p);
%constant double mcmclib_negative_binomial_lpdf (gsl_vector* x, void* in_p);


typedef struct {
	double * p;
	int * n;
} pascal_lpdf_p;

pascal_lpdf_p* mcmclib_pascal_lpdf_alloc (double * p, int * n);
void mcmclib_pascal_lpdf_free (pascal_lpdf_p* p);
%constant double mcmclib_pascal_lpdf (gsl_vector* x, void* in_p);


typedef struct {
	double * p;
} geometric_lpdf_p;

geometric_lpdf_p* mcmclib_geometric_lpdf_alloc (double * p);
void mcmclib_geometric_lpdf_free (geometric_lpdf_p* p);
%constant double mcmclib_geometric_lpdf (gsl_vector* x, void* in_p);


typedef struct {
	double * p;
} logarithmic_lpdf_p;

logarithmic_lpdf_p* mcmclib_logarithmic_lpdf_alloc (double * p);
void mcmclib_logarithmic_lpdf_free (logarithmic_lpdf_p* p);
%constant double mcmclib_logarithmic_lpdf (gsl_vector* x, void* in_p);

%nocallback;
