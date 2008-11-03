dat <- read.csv("ex1_out.csv")
table(dat$proposal)
chain <- dat[,1:2]
chu <- unique(chain)
#acceptance rate:
nrow(chu) / nrow(chain)

nu <- nrow(chu)
plot(chu[1e4 + 1:100,], type="l")
ch2 <- chu[-seq_len(nu/2),]
plot(ch2, type="p", pch=16, cex=0.2)

#load local performance infos
w <- read.csv("ex1_extra_out.csv")
w <- subset(w, proposal != 2)
wx <- w[,1:2]
wp <- w$proposal
wn0 <- 1-rank(w$ntries0[wp==0])/sum(wp==0)
wn1 <- 1-rank(w$ntries1[wp==1])/sum(wp==1)
plot(wn0 ~ wx[wp==0,1], type="n", ylim=c(0, 1))
lines(lowess(wx[wp==0,1], wn0))
abline(h=0.5, lty=2)

w$tag <- apply(w[,c("ntries0", "ntries1")], 1, which.min)
table(w$tag)
plot(x1~x0, col=tag, pch=16, cex=0.2, data=w)
library(MASS)
m.lda <- lda(factor(tag) ~ x0 + x1, data=w)
library(nnet)
m.nnet <- nnet(factor(tag) ~ x0 + x1, data=w, size=2)
set.seed(1234)
x0 <- w[sample(nrow(w), size=5000),c("x0", "x1")]
tag.hat <- as.numeric(predict(m.nnet, x0, type="class"))
plot(x0[,1], x0[,2], col=tag.hat, pch=16, cex=0.2)
plot(x0[,1], x0[,2], col=predict(m.lda, x0)$class, pch=16, cex=0.2)

library(lattice)
xyplot(ntries0~x0 | factor(proposal), data=w)

nrow(w)
w <- subset(w, is.finite(w))
nrow(w)
#library(nnet)
#m0 <- nnet(w ~ x0 + x1, data=w, subset= proposal==0, size=2, linout=TRUE)
#m1 <- nnet(w ~ x0 + x1, data=w, subset= proposal==1, size=2, linout=TRUE)
library(mgcv)
m0 <- gam(w ~ s(x0) + s(x1), data=w, subset= proposal==0)
m1 <- gam(w ~ s(x0) + s(x1), data=w, subset= proposal==1)

y0 <- matrix(, 5000, 2)
y0[, 1] <- predict(m0, newdata=x0)
y0[, 2] <- predict(m1, newdata=x0)
summary(y0)
tags <- apply(y0, 1, which.max)
table(tags)
plot(x0[,1], y0[,1], col=tags, pch=16, cex=0.2)
plot(x0[,2], y0[,1], col=tags, pch=16, cex=0.2)
plot(x0[,1], y0[,2], col=tags, pch=16, cex=0.2)
plot(x0[,2], y0[,2], col=tags, pch=16, cex=0.2)

