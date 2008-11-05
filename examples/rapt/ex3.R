d <- read.csv("ex3_extra_out.csv")
plot(d$score, type="b", xlab="batch number", ylab=expression(theta))
abline(h=0, lty=2)
