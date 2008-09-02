#include "lpdf_hierarchical.h"

post_lpdf_p* mcmclib_lpdf_post_alloc(gsl_vector* x, distrfun_p prior, void* parms,
	distrfun_p loglik, gsl_vector** childs, void** child_parms) {
	post_lpdf_p* ans = (post_lpdf_p*) malloc(sizeof(post_lpdf_p));
	ans->x = x;
	ans->prior = prior;
	ans->parms = parms;
	ans->loglik = loglik;
	ans->childs = childs;
	ans->child_parms = child_parms;
	ans->workspace = gsl_vector_alloc(x->size);
	return ans;
}

void mcmclib_lpdf_post_free(post_lpdf_p* p) {
	gsl_vector_free(p->workspace);
	free(p);
}

double mcmclib_lpdf_post(void* data, gsl_vector* x) {
	post_lpdf_p* d = (post_lpdf_p*) data;
	double ans = 0.0;
	int i=0;
	/*store old value*/
	gsl_vector_memcpy(d->workspace, d->x);
	/*set new value*/
	gsl_vector_memcpy(d->x, x);

	ans += d->prior(d->x, d->parms);
	if(!isfinite(ans)) {
		/*restore old value*/
		gsl_vector_memcpy(d->x, d->workspace);
		return(ans);
	}
	for(int i=0; d->childs[i] != NULL; i++) {
		ans += d->loglik(d->childs[i], d->child_parms[i]);
		if(!isfinite(ans))
			break;
	}

	/*restore old value*/
	gsl_vector_memcpy(d->x, d->workspace);

	return ans;
}
