x <- read.csv("chain.csv")
names(x) <- c("mu", "sigma", paste("k", seq_len(100), sep=""))
library(coda)
xx <- mcmc(x, thin=10)
xx[,2] <- exp(xx[,2])
par(mfrow=c(2,5), mar=c(2,2,2,0))
K <- 10
for(k in seq_len(K)) {
  xk <- xx[,k + 2]
  fk <- table(xk) / length(xk)
  kk <- as.numeric(names(fk))
  plot(kk, fk, type="h", main="", xlab="", ylab="", xlim=c(-8, 8), ylim=c(0, 0.31),
       axes=FALSE)
  axis(1, at=seq(-8, 8))
  axis(2)
}
