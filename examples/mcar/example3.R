P <- 3
DIM <- 95
THIN <- 10

info <- dget("adiacenze.txt")
W <- info$W
num <- info$num

D <- dget("iperischecardcron.txt")
y <- c(D$Y10, D$Y20, D$Y30)
offset <- c(D$N10, D$N20, D$N30)
id <- rep(1:3-1, DIM) * DIM + 1 + rep((1:DIM -1), each=3)
y <- y[id]
offset <- offset[id]

X <- matrix(, DIM*P, P)
for(i in 1:DIM)
  X[(i-1)*P + 1:P, 1:P] <- diag(P)
m <- glm(y~X-1, family=poisson, offset=log(offset))

write.table(y, file="y_2.dat", row.names=FALSE, col.names=FALSE)
write.table(offset, file="offset_2.dat", row.names=FALSE, col.names=FALSE)
write.table(W, file="W_2.dat", row.names=FALSE, col.names=FALSE)

library(coda)
beta <- mcmc(matrix(read.table("chain_beta.dat")[[1]],
                  byrow=TRUE, ncol=P), thin=THIN)
plot(beta)

lpdf <- mcmc(matrix(read.table("chain_lpdf.dat")[[1]],
                    byrow=TRUE, ncol=1), thin=THIN)
plot(lpdf)

g <- function(x) (pi/2) * (exp(x) - 1) / (exp(x) + 1)
as <- mcmc(matrix(read.table("chain_alphasigma.dat")[[1]],
                  byrow=TRUE, ncol=P*(P-1)/2 + P), thin=THIN)
as[,seq_len(P*(P-1)/2)] <- g(as[,seq_len(P*(P-1)/2)])
as[,P*(P-1)/2 + 1:P] <- exp(as[,P*(P-1)/2 + 1:P])
plot(as[,4:6])

a12s <- mcmc(matrix(read.table("chain_alpha12sigma.dat")[[1]],
                  byrow=TRUE, ncol=P*P), thin=THIN)
a12s[,seq_len(P*(P-1))] <- g(a12s[,seq_len(P*(P-1))])
a12s[,P*(P-1) + 1:P] <- exp(a12s[,P*(P-1) + 1:P])
plot(a12s[,7:9])
summary(a12s[,7:9])

phi <- mcmc(matrix(read.table("chain_phi.dat")[[1]],
                  byrow=TRUE, ncol=DIM*P), thin=THIN)
round(apply(phi, 2, var), 3)[1:24]
plot(phi[,1:3])

gg <- mcmc(matrix(read.table("chain_gammaii.dat")[[1]], byrow=TRUE, ncol=P), thin=THIN)
plot(gg)
bb <- mcmc(matrix(read.table("chain_bii.dat")[[1]], byrow=TRUE, ncol=P), thin=THIN)
plot(bb)
