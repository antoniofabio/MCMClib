P <- 3
as <- matrix(read.table("chain_alpha12sigma.dat")[[1]], ncol=P*P)
plot(as[,1], type="l")
plot(as[,P*(P-1)+1], type="l")
plot(as[,P*(P-1)+2], type="l")
plot(as[,P*(P-1)+3], type="l")
