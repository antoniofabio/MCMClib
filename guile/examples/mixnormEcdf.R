library(mvtnorm)
rmixnorm <- function(n, beta, muList, SigmaList) {
  beta <- beta
  K <- length(beta)
  stopifnot(K == length(muList))
  stopifnot(K == length(SigmaList))
  coin <- sample(K, size=n, replace=TRUE, prob=beta)
  ans <- matrix(NA, n, length(muList[[1]]))
  for(i in seq_len(n))
    ans[i,] <- rmvnorm(1, mean=muList[[coin[i]]], sigma=SigmaList[[coin[i]]])
  return(ans)
}
N <- 1000
set.seed(1234)
X0 <- rmixnorm(N, beta=rep(1, 2),
               muList=list(rep(-1, 5), rep(1, 5)),
               SigmaList=list(diag(1, 5, 5), diag(1, 5, 5)))
write.table(X0, file="mixnorm_iid_sample.dat", sep=" ", dec=".",
            col.names=FALSE, row.names=FALSE)
