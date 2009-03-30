N <- 9
MU <- 1.5
RHO <- 2.0
SIGMA2 <- 2.0
TAUSQ <- 0.5

xy <- matrix(, N, 2)
k <- 1
for(i in 0:2) for(j in 0:2) {
  xy[k, 1] <- i
  xy[k, 2] <- j
  k <- k + 1
}
D <- as.matrix(dist(xy))

library(mvtnorm)
getSigma <- function(rho, sigma, tausq, D) {
  N <- nrow(D)
  Sigma <- diag(sigma, N)
  for(i in 1:(N-1)) for(j in (i+1):N) {
    dij <- D[i,j]
    cij <- (SIGMA2 - TAUSQ) * ( 1 - exp(- dij / RHO) )
    if(dij > 0)
      cij <- cij + TAUSQ
    Sigma[i,j] <- Sigma[j,i] <- SIGMA2 - cij
  }
  return(Sigma)
}
Sigma <- getSigma(RHO, SIGMA2, TAUSQ, D)

set.seed(1234)
Y <- rmvnorm(1, rep(0, N), Sigma)

##read chain output
X <- read.csv("chain.csv")
library(coda)
xx <- mcmc(exp(X[,1:3]))
plot(xx)
