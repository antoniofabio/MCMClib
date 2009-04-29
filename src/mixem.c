/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#include<assert.h>
#include<gsl/gsl_vector.h>
#include<gsl/gsl_matrix.h>
#include "mixem.h"
#include "mvnorm.h"
#include "mixnorm.h"
#include "mixem_rec.h"

/**Fitting a gaussian mixture distribution by the EM algorithm
*/
void mcmclib_mixem_fit(gsl_matrix* X,
		       gsl_vector** mu, gsl_matrix** Sigma,
		       gsl_vector* w, int NITER) {
  int K = w->size;
  int N = X->size1;

  for(int ITER=0; ITER < NITER; ITER++) {
    mcmclib_mixem_rec* m = mcmclib_mixem_rec_alloc(mu, Sigma, w);

    /*accumulate rows info in object 'm'*/
    for(int n=0; n<N; n++) {
      gsl_vector_view rv = gsl_matrix_row(X, n);
      mcmclib_mixem_rec_add(m, &(rv.vector));
    }
    
    /*update means and variances estimates*/
    mcmclib_mixem_rec_update(m);

    for(int k=0; k<K; k++)
      assert(gsl_vector_get(w, k) <= 1.0);

    mcmclib_mixem_rec_free(m);
  }
}
