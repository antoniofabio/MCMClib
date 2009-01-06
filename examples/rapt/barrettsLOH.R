X <- read.csv("barrettsLOH_X.csv")
X <- as.matrix(X)
N <- 50000
K <- 20
DIM <- 4
XX <- array(X, c(DIM, K, N+1))
dimnames(XX) <- list(par=colnames(X), chain=1:K, time=NULL)
X1 <- t(XX[,1,])
plot(as.ts(X1[1:35000,]))
library(coda)

