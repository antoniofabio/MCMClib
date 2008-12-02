#include <stdio.h>
#include <gsl/gsl_rng.h>
#include <mixolem_suff.h>
#include "mixnorm_target.h"

/**fitted mixture parmaters*/
extern mcmclib_mixolem_suff* gamma_hat;
/*output files handlers*/
extern FILE *out_beta, *out_mu, *out_Sigma;

/*init fitted mixture data*/
void mixnorm_fitting_init(gsl_rng* r);
void mixnorm_fitting_free();
void mixnorm_fitting_store();
