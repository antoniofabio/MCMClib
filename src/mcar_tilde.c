#include <assert.h>
#include <math.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_blas.h>
#include "mcar_tilde.h"

mcmclib_mcar_tilde_lpdf* mcmclib_mcar_tilde_lpdf_alloc(int p, int n, gsl_matrix* M) {
	mcmclib_mcar_tilde_lpdf* a = (mcmclib_mcar_tilde_lpdf*) malloc(sizeof(mcmclib_mcar_tilde_lpdf));
	assert(p>0);
	assert(n>0);
	a->p = p;
	a->n = n;
	a->B_tilde = gsl_matrix_alloc(p, p);
	a->alpha1 = gsl_vector_alloc(p * (p-1) / 2);
	a->alpha2 = gsl_vector_alloc(p * (p-1) / 2);
	a->sigma = gsl_vector_alloc(p);
	gsl_vector_set_all(a->alpha1, 0.0);
	gsl_vector_set_all(a->alpha2, 0.0);
	gsl_vector_set_all(a->sigma, 0.4);
	gsl_vector_set(a->sigma, 0, 0.5);
	a->M = gsl_matrix_alloc(n, n);
	gsl_matrix_memcpy(a->M, M);
	a->m = gsl_vector_alloc(n);
	for(int i=0; i<n; i++) {
		int count = 0;
		for(int j=0; j<n; j++)
			count += gsl_matrix_get(M, i, j) == 1.0;
		gsl_vector_set(a->m, i, (double) count);
	}

	a->vcov = gsl_matrix_alloc(p*n, p*n);
	gsl_matrix_set_identity(a->vcov);
	a->mu = gsl_vector_alloc(p*n);
	gsl_vector_set_zero(a->mu);
	a->mvnorm = mcmclib_mvnorm_lpdf_alloc(a->mu, a->vcov->data);
	return a;
}

void mcmclib_mcar_tilde_lpdf_free(mcmclib_mcar_tilde_lpdf* p) {
	mcmclib_mvnorm_lpdf_free(p->mvnorm);
	gsl_matrix_free(p->vcov);
	gsl_vector_free(p->mu);
	gsl_vector_free(p->m);
	gsl_matrix_free(p->M);
	gsl_vector_free(p->sigma);
	gsl_vector_free(p->alpha1);
	gsl_vector_free(p->alpha2);
	gsl_matrix_free(p->B_tilde);
  free(p);
}

void mcmclib_mcar_tilde_lpdf_set_alpha(mcmclib_mcar_tilde_lpdf* p,
																			 gsl_vector* alpha1, gsl_vector* alpha2) {
	gsl_vector_memcpy(p->alpha1, alpha1);
	gsl_vector_memcpy(p->alpha2, alpha2);
}

void mcmclib_mcar_tilde_lpdf_set_sigma(mcmclib_mcar_tilde_lpdf* p,
																			 gsl_vector* sigma) {
	gsl_vector_memcpy(p->sigma, sigma);
}

static int ALPHA(int i, int j) {
	assert(j > i);
	int ans = 1;
	for(int i1 = 0; i1<i; i1++)
		for(int j1 = i1+1; j1<j; j1++)
			ans++;
	return ans;
}

static double alpha_get(gsl_vector* a, int i, int j) {
	return gsl_vector_get(a, ALPHA(i, j));
}

static void Givens_set_Shij(gsl_matrix* S, int i, int j, double alpha_ij) {
	gsl_matrix_set_identity(S);
	gsl_matrix_set(S, i, i, cos(alpha_ij));
	gsl_matrix_set(S, j, j, cos(alpha_ij));
	gsl_matrix_set(S, i, j, sin(alpha_ij));
	gsl_matrix_set(S, j, i, -sin(alpha_ij));
}

static void Givens_rotations(gsl_matrix* A, gsl_vector* alpha) {
	assert(A->size1 == A->size2);
	int p = A->size1;
	assert(alpha->size == (p * (p-1) / 2));
	gsl_matrix* S[2];
	for(int h=0; h<2; h++)
		S[h] = gsl_matrix_alloc(p, p);
	gsl_matrix_set_zero(A);
	gsl_matrix_set_identity(S[1]);
	
	for(int i=0; i<(p-1); i++) {
		for(int j=i+1; j<p; j++) {
			Givens_set_Shij(S[0], i, j, alpha_get(alpha, i, j));
			gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, S[1], S[0], 0.0, A);
			gsl_matrix_memcpy(S[1], S[2]);
		}
	}
	
	for(int h=0; h<2; h++)
		gsl_matrix_free(S[h]);
}

/* rebuild asymm. matrix from its SVD decomposition */
static void anti_SVD(gsl_matrix* A, gsl_matrix* P1, gsl_matrix* P2, gsl_vector* sigma) {
	int p = A->size1;
	gsl_matrix* Delta = gsl_matrix_alloc(p, p);
	gsl_matrix* P1Delta = gsl_matrix_alloc(p, p);
	gsl_matrix_set_zero(Delta);
	for(int i=0; i<p; i++)
		gsl_matrix_set(Delta, i, i, gsl_vector_get(sigma, i));
	gsl_matrix_set_zero(P1Delta);
	gsl_blas_dgemm(CblasNoTrans, CblasNoTrans, 1.0, P1, Delta, 0.0, P1Delta);
	gsl_blas_dgemm(CblasNoTrans, CblasTrans, 1.0, P1Delta, P2, 0.0, A);
	gsl_matrix_free(Delta);
	gsl_matrix_free(P1Delta);
}

static void get_B_tilde(gsl_matrix* A, gsl_vector* sigma,
												gsl_vector* alpha1, gsl_vector* alpha2) {
	int p = A->size1;
	gsl_matrix* P1 = gsl_matrix_alloc(p, p);
	gsl_matrix* P2 = gsl_matrix_alloc(p, p);
	Givens_rotations(P1, alpha1);
	Givens_rotations(P2, alpha2);
	anti_SVD(A, P1, P2, sigma);
	gsl_matrix_free(P1);
	gsl_matrix_free(P2);
}

static int is_positive_definite(mcmclib_mcar_tilde_lpdf* p) {
	double sigma_0 = gsl_vector_get(p->sigma, 0);
	if((sigma_0 <= 0.0) ||(sigma_0 >= 1.0))
		return 0;

	get_B_tilde(p->B_tilde, p->sigma, p->alpha1, p->alpha2);
	for(int i=0; i<p->n; i++) {
		double mi = gsl_vector_get(p->m, i);
		double mUi = 0.0;
		for(int j=0; j<p->n; j++)
			if(i < j)
				mUi += gsl_matrix_get(p->M, i, j);
		double mLi = mi - mUi;
		
		for(int k=0; k<p->p; k++) {
			double a = fabs(gsl_matrix_get(p->B_tilde, k, k));
			double b = 0.0;
			for(int l=0; l<p->p; l++) {
				if(l == k)
					continue;
				b += fabs(gsl_matrix_get(p->B_tilde, k, l));
			}
			b *= mUi / mi;
			double c = 0.0;
			for(int l=0; l<p->p; l++) {
				if(l == k)
					continue;
				c += fabs(gsl_matrix_get(p->B_tilde, l, k));
			}
			c *= mLi / mi;

			if((a + b + c >= 1.0))
				return 0;
		}
	}
	return 1;
}

double mcmclib_mcar_tilde_lpdf_compute(void* in_p, gsl_vector* x) {
  mcmclib_mcar_tilde_lpdf* p = (mcmclib_mcar_tilde_lpdf*) in_p;
	if(!is_positive_definite(p))
		return log(0.0);
	gsl_matrix_set_identity(p->vcov); /*FIXME*/
  return mcmclib_mvnorm_lpdf_compute(p->mvnorm, x);
}
