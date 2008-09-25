fn <- as.list(paste("chain_", 0:19,".csv", sep=""))
dd <- lapply(fn, read.csv)
save.image("LOH_chains_20.RData")
load("LOH_chains_20.RData")
library(coda)
chains <- do.call(mcmc.list, lapply(dd, mcmc))
rejectionRate(chains)
summary(chains)

chain <- do.call(rbind, chains)
set.seed(1234)
pi <- chain[sample(nrow(chain), 5e4, replace=TRUE), 2:3]
pdf("example3_pi1_pi2.pdf")
plot(pi, pch=16, cex=0.5, col=rgb(0,0,0, 0.05),
	xlab=expression(pi[1]), ylab=expression(pi[2]))
dev.off()
