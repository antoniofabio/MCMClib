source("mcmclib.R")
attach(rawgsl)
attach(mcmclib)

dyn.load("/usr/lib/libgslcblas.so", local=FALSE)
dyn.load("/usr/lib/libgsl.so", local=FALSE)
x <- rnorm(23)
gx <- vec2gvec(x)
x1 <- gvec2vec(gx)
stopifnot(identical(x, x1))
gsl_vector_free(gx)

pp <- getNativeValue(getNativeSymbolInfo("gsl_rng_default")$address,
                     pointerType)
rng <- gsl_rng_alloc(pp, RETURN=pointerType)

xx <- replicate(1000, gsl_ran_gaussian(rng, 1.0, RETURN=doubleType))

dyn.load("../src/libmcmclib.so", local=FALSE)
dim <- 1000
x <- vec2gvec(rep(0.0, dim))

dyn.load("f.so")
f <- getNativeSymbolInfo("f")$address
sampler <- mcmclib_gauss_scalar_am_alloc(rng, f, 0L, x, 0.1, 1000L,
                                         RETURN = pointerType)
system.time(mcmclib_amh_update_N(sampler, 10000L))
y <- gvec2vec(x)
mcmclib_amh_free(sampler)

pp <- getNativeValue(getNativeSymbolInfo("gsl_rng_default")$address,
                     pointerType)
rng <- gsl_rng_alloc(pp, RETURN=pointerType)
gvec_set(x, seq_len(dim)-1, 0.0)
sampler <- mcmclib_gauss_rw_alloc(rng, f, 0L, x, 0.1,
                                  RETURN = pointerType)
system.time(mcmclib_mh_update_N(sampler, 10000L))
y <- gvec_get(x, seq_len(dim)-1)
