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
