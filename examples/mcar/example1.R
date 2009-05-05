P <- 3
THIN <- 10
as <- ts(matrix(read.table("chain_alpha12sigma.dat")[[1]],
                byrow=TRUE, ncol=P*P), freq=1/THIN)
ls <- ts(read.table("chain_lpdf.dat")[[1]], freq=1/THIN)
plot(ls)
plot(as[,1])
offset <- P*(P-1)
sigma <- as[,offset+1:3]
plot(sigma)
pairs(sigma, pch=16)
