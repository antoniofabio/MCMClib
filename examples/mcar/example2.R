P <- 3
THIN <- 10
g <- function(x) (pi/2) * (exp(x) - 1) / (exp(x) + 1)

library(coda)
beta <- mcmc(matrix(read.table("chain_beta.dat")[[1]],
                    byrow=TRUE, ncol=P), thin=THIN)
plot(beta)
