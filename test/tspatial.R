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
lpdf <- function(s) {
  Sigma <- diag(SIGMA2, N)
  for(i in 1:(N-1)) for(j in (i+1):N) {
    dij <- D[i,j]
    cij <- (SIGMA2 - TAUSQ) * ( 1 - exp(- dij / RHO) )
    if(dij > 0)
      cij <- cij + TAUSQ
    Sigma[i,j] <- Sigma[j,i] <- SIGMA2 - cij
  }
  print(eigen(Sigma)$values)
  dmvnorm(rep(s, N), rep(MU, N), Sigma, log=TRUE)
}
lpdf(1)
lpdf(0)
lpdf(-1)
lpdf(1.5)
