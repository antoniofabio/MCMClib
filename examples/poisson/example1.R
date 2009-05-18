xx <- matrix(read.table("chain_beta.dat")[[1]], ncol=3, byrow=TRUE)

library(coda)
XX <- mcmc(xx, thin=10)
plot(XX)

X <- rbind(diag(3), matrix(0, 3, 3))
y <- rep(3, 6)
m <- glm(y ~ X - 1, family=poisson)
summary(m)
