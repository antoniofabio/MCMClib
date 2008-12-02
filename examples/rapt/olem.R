DIM <- 5
N <- 25000
M <- 4
X <- array(read.table("X.dat")[[1]], dim=c(DIM, M, N))
X1 <- X[1,,]
plot(as.ts(t(X1)))

mu <- list(rep(-1.5, DIM), rep(1.5, DIM))
Sigma <- list(matrix(-0.1, DIM, DIM), matrix(-0.1, DIM, DIM))
diag(Sigma[[1]]) <- diag(Sigma[[2]]) <- 1
Sigma[[1]] <- Sigma[[1]] * 4
beta <- list(0.2, 0.8)
rpi <- function(n) {
  require(mvtnorm) || stop()
  N1 <- round(n * beta[[1]])
  N2 <- n - N1
  Y1 <- rmvnorm(N1, mu[[1]], Sigma[[1]])
  Y2 <- rmvnorm(N2, mu[[2]], Sigma[[2]])
  Y <- rbind(Y1, Y2)
  Y[sample(n),]
}
xx <- rpi(5000)[,1]
pi <- density(xx, from=-8, to=6)
par(mar=c(2,4,2,2))
hist(X[1,,], breaks=100, prob=TRUE, col="gray", xlab="", main="")
lines(pi, lwd=2, lty=2)
