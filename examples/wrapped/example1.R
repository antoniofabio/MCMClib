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
  lL <- c(-8, 8)
  plot(kk, fk, type="h", main="", xlab="", ylab="", xlim=lL, ylim=c(0, 0.70),
       axes=FALSE)
  axis(1, at=seq(lL[1], lL[2]))
  axis(2)
  box()
}
