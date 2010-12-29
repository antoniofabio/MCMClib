source("mcmclib.R")
attach(rawgsl)
attach(mcmclib)

dyn.load("/usr/lib/libgslcblas.so", local=FALSE)
dyn.load("/usr/lib/libgsl.so", local=FALSE)
dyn.load("../src/libmcmclib.so", local=FALSE)

cv <- function(x, y) mean(x*y) - mean(x)*mean(y)
v2 <- function(x) cv(x, x)

g1 <- function(xx, lag.max) {
  xx <- as.matrix(xx)
  m <- mcmclib_monitor_acf_alloc(NCOL(xx), as.integer(lag.max),
                                 RETURN = pointerType)
  for(i in seq_len(NROW(xx))) {
    mcmclib_monitor_acf_update(m, vec2gvec(xx[i,]))
  }
  acf <- mat2gmat(matrix(0.0, lag.max+1, NCOL(xx)))
  mcmclib_monitor_acf_get(m, acf)
  mcmclib_monitor_acf_free(m)
  ans <- gmat2mat(acf)
  gsl_matrix_free(acf)
  return(ans)
}

g <- function(x, ...) acf(x, plot=FALSE, ...)$acf

ntail <- function(x, n) if(n==0) x else tail(x, -n)
nhead <- function(x, n) if(n==0) x else head(x, -n)

g2 <- function(x, l) {
  n <- length(x)
  m <- mean(x)
  x.l <- nhead(x, l) ## last n-l obs
  x.sl <- ntail(x, l) ## first n-l obs
  s.xy <- sum(x.l * x.sl)
  s.l <- sum(x.l)
  s.sl <- sum(x.sl)
  C <- (s.xy - m * ( s.l + s.sl ) + (n-l) * m^2)/n
  message("covariance = ", round(C, 3))
  return(C/v2(x))
}

set.seed(1234)
xx <- rnorm(10)
g1(xx, lag.max=0)

set.seed(1234)
xx <- arima.sim(n = 1000,
                list(ar = c(-0.3, 0.3), ma = c()),
                sd = sqrt(2))
g(xx, lag.max=3)
g1(xx, lag.max=3)

##
## Now test the multivariate case
##
set.seed(1234)
n <- 1000
gen <- function() arima.sim(n=n, list(ar=c(-0.3, 0.3), ma=c()), sd = sqrt(2))
xx <- cbind(x1=gen(), x2=gen(), x3=gen())

g(xx, lag.max=3)
g1(xx, lag.max=3)

##
## Test IACT
##
set.seed(1234)
n <- 1000
gen <- function() arima.sim(n=n, list(ar=0.95, ma=c()), sd = sqrt(2))
xx <- cbind(x1=gen(), x2=gen(), x3=gen())

A <- g1(xx, lag.max=100)
iact_from_acf(A)
