/*
 *  MCMClib: A C Library for doing MCMC
 *  Copyright (C) 2009 Antonio, Fabio Di Narzo
 *
 *  This program is free software; you can redistribute it and/or modify
 *  it under the terms of the GNU General Public License as published by
 *  the Free Software Foundation; either version 2 of the License, or
 *  (at your option) any later version.
 */
#ifndef __INCA_H__
#define __INCA_H__

/**\addtogroup adaptive
@{
\defgroup inca INter-Chain Adaptation
@{*/

#include "amh.h"

/**\brief Generic INCA sampler */
typedef struct {
  mcmclib_amh* amh;
  gsl_vector** x; /**< array of current chain values*/
  int M; /**< number of parallel chains*/
} mcmclib_inca;

/**\brief alloc a new INCA sampler object
@param amh \ref mcmclib_amh object, allocated with \ref mcmclib_amh_alloc
@param x array of M current chain values
@param M number of parallel chains
*/
mcmclib_inca* mcmclib_inca_alloc(mcmclib_amh* amh, gsl_vector** x, int M);

/**\brief Free a previously allocated INCA sampler object*/
void mcmclib_inca_free(mcmclib_inca* p);

/**\brief update all chains values sequentially
\return 1 if last chain move accepted, 0 if rejected*/
int mcmclib_inca_update(mcmclib_inca* p);

/**@}*/
/**@}*/

#endif
