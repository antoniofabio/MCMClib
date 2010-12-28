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
                 pst("rng_", c("alloc", "free", "set")),
                 pst("ran_", pst(distribs, c("", "_pdf"))))),
           rawgsl)

dyn.constant(pst("gsl_rng_",
                 c("default",
                   "mt19937", "ranlxs0", "ranlxs1", "ranlxs2",
                   "ranlxd1", "ranlxd2", "ranlux", "ranlux389",
                   "cmrg", "mrg", "taus", "taus2", "gfsr4")),
             rawgsl)

## utility functions and additional types
with(rawgsl, {

  gvec.type <- structType(list(size = uint64Type, stride = uint64Type,
                               data = pointerType, block = pointerType,
                               owner = sint32Type))

  gvec_size <- function(x) {
    getStructField(x, "size", gvec.type)
  }
  
  gvec_get <- function(x, ii) sapply(as.list(ii), function(i) gsl_vector_get(x, as.integer(i), RETURN=doubleType))

  gvec2vec <- function(x) gvec_get(x, seq_len(gvec_size(x))-1)

  gvec_set <- function(x, ii, val) {
    mapply(function(i, val.i) gsl_vector_set(x, as.integer(i), as.double(val.i)),
           as.list(ii), as.list(val))
    return(invisible(val))
  }

  vec2gvec <- function(x) {
    y <- gsl_vector_alloc(as.integer(length(x)), RETURN = pointerType)
    gvec_set(y, seq_len(length(x))-1, x)
    return(y)
  }
  
  gmat.type <- structType(list(size1 = uint64Type, size2 = uint64Type,
                               tda = uint64Type,
                               data = pointerType, block = pointerType,
                               owner = sint32Type))

  gmat_size1 <- function(x) {
    getStructField(x, "size1", gmat.type)
  }
  gmat_size2 <- function(x) {
    getStructField(x, "size2", gmat.type)
  }
  gmat_dim <- function(x) {
    c(gmat_size1(x), gmat_size2(x))
  }

  gmat_get <- function(x, ii, jj) {
    outer(ii, jj, Vectorize(function(i, j) {
      gsl_matrix_get(x, as.integer(i), as.integer(j), RETURN=doubleType)
    }))
  }

  gmat2mat <- function(x) {
    gmat_get(x, seq_len(gmat_size1(x))-1, seq_len(gmat_size2(x))-1)
  }

  gmat_set <- function(x, i, j, value) {
    gsl_matrix_set(x, as.integer(i), as.integer(j), as.double(value))
  }
  
  mat2gmat <- function(x) {
    y <- gsl_matrix_alloc(0.0, NROW(x), NCOL(x), RETURN = pointerType)
    for(i in seq_len(NROW(x))) for(j in seq_len(NCOL(x))) {
      gmat_set(y, i-1, j-1, x[i, j])
    }
    return(y)
  }
  
})
