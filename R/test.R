source("mcmclib.R")
attach(rawgsl)
attach(mcmclib)

dyn.load("/usr/lib/libgslcblas.so", local=FALSE)
dyn.load("/usr/lib/libgsl.so", local=FALSE)
x <- gsl_vector_alloc(23L, RETURN=pointerType)
gsl_vector_set(x, 7L, -2)
gsl_vector_get(x, 7L, RETURN=doubleType)
gsl_vector_free(x)

pp <- getNativeValue(getNativeSymbolInfo("gsl_rng_default")$address,
                     pointerType)
rng <- gsl_rng_alloc(pp, RETURN=pointerType)

xx <- replicate(1000, gsl_ran_gaussian(rng, 1.0, RETURN=doubleType))

dyn.load("../src/libmcmclib.so", local=FALSE)
dim <- 1000L
x <- gsl_vector_alloc(dim, RETURN=pointerType)
gvec_set(x, seq_len(dim)-1, 0.0)

dyn.load("f.so")
f <- getNativeSymbolInfo("f")$address
sampler <- mcmclib_gauss_scalar_am_alloc(rng, f, 0L, x, 0.1, 1000L,
                                         RETURN = pointerType)
system.time(mcmclib_amh_update_N(sampler, 100000L))
uh <- gvec_get(x, seq_len(dim)-1)
mcmclib_amh_free(sampler)

pp <- getNativeValue(getNativeSymbolInfo("gsl_rng_default")$address,
                     pointerType)
rng <- gsl_rng_alloc(pp, RETURN=pointerType)
gvec_set(x, seq_len(dim)-1, 0.0)
sampler <- mcmclib_gauss_rw_alloc(rng, f, 0L, x, 0.1,
                                  RETURN = pointerType)
system.time(mcmclib_mh_update_N(sampler, 100000L))
uh <- gvec_get(x, seq_len(dim)-1)
