d <- read.csv("ex3_extra_out.csv")
plot(d$score, type="b", xlab="batch number", ylab=expression(theta))
abline(h=0, lty=2)

ths <- seq(-3.0, 3.0, by=0.5)
score <- numeric()
system.time(
for(th in ths) {
  message("evaluating th = ", th)
  args <- paste(th, "100000 2")
  score[as.character(th)] <- as.numeric(system(paste("./ex3", args), intern=TRUE))
})

pdf("tmp.pdf")
plot(ths, score, type="b", xlab="threshold", ylab="acceptance rate")
dev.off()
