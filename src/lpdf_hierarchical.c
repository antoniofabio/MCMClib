#include "lpdf_hierarchical.h"

post_lpdf_p* mcmclib_lpdf_post_alloc(distrfun_p prior, void* parms,
	distrfun_p loglik, gsl_vector** childs, void** child_parms) {
	post_lpdf_p* ans = (post_lpdf_p*) malloc(sizeof(post_lpdf_p));
	ans->prior = prior;
	ans->parms = parms;
	ans->loglik = loglik;
	ans->childs = childs;
	ans->child_parms = child_parms;
	return ans;
}

void mcmclib_lpdf_post_free(post_lpdf_p* p) {
	free(p);
}

double mcmclib_lpdf_post(gsl_vector* x, void* data) {
	post_lpdf_p* d = (post_lpdf_p*) data;
	double ans = 0.0;
	int i=0;

	ans += d->prior(x, d->parms);
	for(int i=0; d->childs[i] != NULL; i++)
		ans += d->loglik(d->childs[i], d->child_parms[i]);

	return ans;
}
