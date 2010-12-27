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
x <- gsl_vector_alloc(100L, RETURN=pointerType)
gvec_set(x, 0:99, 0.0)

dyn.load("f.so")
f <- getNativeSymbolInfo("f")$address
sampler <- mcmclib_gauss_scalar_am_alloc(rng, f, 0L, x, 3.0, 1000L,
                                         RETURN = pointerType)
for(i in seq_len(2e4)) mcmclib_amh_update(sampler)
uh <- Vectorize(function(i) gsl_vector_get(x, as.integer(i), RETURN=doubleType))(seq.int(0, 99))
mcmclib_amh_free(sampler)
