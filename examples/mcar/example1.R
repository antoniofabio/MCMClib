P <- 6
THIN <- 10
offset <- P*(P-1)
offset2 <- P*(P-1)/2
g <- function(x) (pi/2) * (exp(x) - 1) / (exp(x) + 1)

a12s <- ts(matrix(read.table("chain_alpha12sigma.dat")[[1]],
                byrow=TRUE, ncol=P*P), freq=1/THIN)
as <- ts(matrix(read.table("chain_alphasigma.dat")[[1]],
                byrow=TRUE, ncol=offset2 + P), freq=1/THIN)
a12s[,1:offset] <- g(a12s[,1:offset])
sigma <- exp(a12s[,offset+1:P])
as[,1:offset2] <- g(as[,1:offset2])
sigmag <- exp(as[,offset2+1:P])
ls <- ts(read.table("chain_lpdf.dat")[[1]], freq=1/THIN)

plot(ls)
plot(a12s[,1:P])
plot(sigma)
pairs(sigma, pch=16)

plot(as[,1:P])
plot(sigmag)

