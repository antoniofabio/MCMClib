source("lazyffi.R")

pst <- function(x, y) {
  ans <- vector(length(x) * length(y), mode="character")
  for(i in seq_along(x)) for(j in seq_along(y)) {
    ans[(i-1)*length(y) + j] <- paste(x[i], y[j], sep="")
  }
  return(ans)
}

rawgsl <- new.env()

distribs <- c("gaussian", "exponential", "laplace", "cauchy",
              "rayleigh", "exppow", "flat", "lognormal", "chisq",
              "landau", "fdist", "tdist", "beta", "logistic",
              "pareto", "weibull")

dyn.import(pst("gsl_",
               c(pst(c("vector_", "matrix_"), c("alloc", "free", "set", "get")),
                 pst("rng_", c("alloc", "free")),
                 pst("ran_", pst(distribs, c("", "_pdf"))))),
           rawgsl)

## utility functions
with(rawgsl, {
  gvec_get <- function(x, ii) sapply(as.list(ii), function(i) gsl_vector_get(x, as.integer(i), RETURN=doubleType))
  gvec_set <- function(x, ii, val) {
    mapply(function(i, val.i) gsl_vector_set(x, as.integer(i), as.double(val.i)),
           as.list(ii), as.list(val))
    return(invisible(val))
  }
  gmat_get <- function(x, ii, jj) {
    outer(ii, jj, function(i, j) {
      gsl_matrix_get(x, as.integer(i), as.integer(j), RETURN=doubleType)
    })
  }
})
