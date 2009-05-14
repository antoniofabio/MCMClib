/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __MCMCLIB_LM_H__
#define __MCMCLIB_LM_H__

#include "mvnorm.h"
/**\addtogroup Models
 @{*/

/**\addtogroup LinearModel
 @{*/

/**\brief Linear regression model with conjugte priors */
typedef struct {
} mcmclib_lm;

/** Alloc extra data for an lm model */
mcmclib_lm* mcmclib_lm_alloc();

/** Free extra data for an lm object */
void mcmclib_lm_free(mcmclib_lm* p);

/**@}*/
/**@}*/
#endif
