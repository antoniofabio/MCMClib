#ifndef __MIXNORM_TGT__
#define __MIXNORM_TGT__

/*state space dimension*/
extern const int DIM;
/*absolute local mean value*/
extern const double MU0;
/*local variances multipliers*/
extern const double V0[];
/*local pairwise correlations*/
extern const double RHO[];
/*mixture proportion of component 1*/
extern const double BETA;
/*number of mixture components*/
#define K 2

mcmclib_mixnorm_lpdf* mixnorm_target_alloc();
void mixnorm_target_free(mcmclib_mixnorm_lpdf* p);

#endif
