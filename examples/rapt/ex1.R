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
library(lattice)
xyplot(ntries0~x0 | factor(proposal), data=w)

w <- read.csv("ex1_extra_out.csv")
wx <- w[,1:2]
w$n <- with(w, ifelse(proposal == 0, ntries0, ntries1))
w$w <- with(w, jump / n)
coplot(w ~ x0 | factor(proposal), data=w)

coplot(wd ~ wn | factor(proposal), data=w)

##doesn't work!
m <- list()
library(nnet)
m[[1]] <- nnet(w ~ x0 + x1, data=w, subset= proposal==0, size=4, linout=TRUE, maxit=1e3)
m[[2]] <- nnet(w ~ x0 + x1, data=w, subset= proposal==1, size=4, linout=TRUE, maxit=1e3)
y <- lapply(m, predict, newdata=w)
w$y <- do.call("-", y)
plot(y ~ x0, col=(y<0)+1, data=w)
###

wp <- w$proposal
wn0 <- w$ntries0[wp==0]
wn0 <- 1-rank(wn0)/length(wn0)
wn1 <- 1-rank(w$ntries1[wp==1])/sum(wp==1)
plot(wn0 ~ wx[wp==0,1], type="n", ylim=c(0, 1))
lines(lowess(wx[wp==0,1], wn0))
abline(h=0.5, lty=2)

w$w <- ifelse(wp==0, wn0, wn1)
m <- list()
library(mgcv)
m[[1]] <- gam(w ~ s(x0) + s(x1), data=w, subset= proposal==0)
m[[2]] <- gam(w ~ s(x0) + s(x1), data=w, subset= proposal==1)
w$w0 <- predict(m[[1]], newdata=w)
w$w1 <- predict(m[[2]], newdata=w)

w$tag <- apply(w[,c("w0", "w1")], 1, which.min)
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
#uh, unexpected result
