chain <- as.ts(read.csv("data.csv"))

plot(window(chain, end=300), main="")
abline(v=1000, lty=2)

pdf("AM_0-100.pdf")
acf(window(chain, end=100))
dev.off()
pdf("AM_101-200.pdf")
acf(window(chain, 101, 200))
dev.off()
pdf("AM_201-300.pdf")
acf(window(chain, 201, 300))
dev.off()

