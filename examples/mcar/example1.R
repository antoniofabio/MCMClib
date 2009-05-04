P <- 3
as <- matrix(read.table("chain_alpha12sigma.dat")[[1]], ncol=P*P)
as <- as[seq(1, 5000, by=100),]
plot(as[,1], type="l")
plot(as[,P*(P-1)+1], type="l")
