d <- read.csv("ex2_extra_out.csv")
plot(d$th, type="b", xlab="batch number", ylab=expression(theta), ylim=c(-1, 1))
abline(h=0, lty=2)
