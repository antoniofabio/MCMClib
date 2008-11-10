X <- matrix(read.table("tmixmem_X.dat")$V1,, 2, byrow=TRUE)
hist(X[,1], breaks=60)
hist(X[,2], breaks=60)

