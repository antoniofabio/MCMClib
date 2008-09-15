filenames <- as.list(paste("data_", 0:2, ".csv", sep=""))

library(coda)
chains <- lapply(filenames, function(i) mcmc(read.csv(i)))
chains <- do.call(mcmc.list, chains)

pdf("INCA_1-100.pdf")
plot(window(chains, end=100), main="")
dev.off()
pdf("INCA_99900-1e5.pdf")
plot(window(chains, start=99900, end=1e5), main="")
dev.off()

pdf("INCA_density_half.pdf")
densityplot(window(chains, start=50001))
dev.off()

sink("INCA_summary.txt")
summary(chains)
sink()
