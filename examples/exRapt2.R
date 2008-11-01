dat <- read.csv("data_rapt2.csv")
table(dat$proposal)
chain <- dat[,1:2]
chu <- unique(chain)
#acceptance rate:
nrow(chu) / nrow(chain)

nu <- nrow(chu)
plot(chu[1e4 + 1:100,], type="l")
ch2 <- chu[-seq_len(nu/2),]
plot(ch2, type="p", pch=16, cex=0.2)
