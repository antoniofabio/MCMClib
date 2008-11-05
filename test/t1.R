library(mixtools)

DIM <- 10
MU0 <- 1.5
V0 <- 2.0

mu <- rep(MU0, DIM)
Sigma <- diag(DIM) * V0 + matrix(1, DIM, DIM) - diag(DIM)

xi <- seq(-3, 3, length=20)
ans <- numeric()
for(i in seq_along(xi))
  ans[i] <- logdmvnorm(rep(xi[i], DIM), mu, Sigma)

write.table(cbind(xi, ans), file="t1.check.dat",
            col.names=FALSE, row.names=FALSE, sep=" ")
