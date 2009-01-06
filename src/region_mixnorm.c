#include "mixnorm.h"
#include "region_mixnorm.h"

int mcmclib_region_mixnorm_compute(gsl_vector* x, void* in_p) {
  mcmclib_mixnorm_lpdf* p = (mcmclib_mixnorm_lpdf*) in_p;
  double pik;
  int ans = 0;
  double pimax = log(0.0);
  for(int k=0; k < p->w->size; k++) {
    pik = mcmclib_mvnorm_lpdf_compute(p->pis[k], x);
    if(pik > pimax) {
      pimax = pik;
      ans = k;
    }
  }
  return ans;
}