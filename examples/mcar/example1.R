P <- 3
DIM <- 95
THIN <- 10

B <- diag(3:1/4)
set.seed(1234)
Gamma <- crossprod(matrix(rnorm(P*P, sd=0.1), P))
W <- matrix(0, DIM, DIM)
for(i in 1:(DIM-1))
  W[i,i+1] <- W[i+1,i] <- 1

asMatrix <- function(x) {
  for(i in seq_len(nrow(x))) {
    row <- x[i, 1][[1]]
    for(j in 1 + seq_len(ncol(x) - 1))
      row <- cbind(row, x[i, j][[1]])
    if(i == 1)
      ans <- row
    else
      ans <- rbind(ans, row)
  }
  return(ans)
}

makePrecision <- function(W, Gamma, B) {
  A <- t(chol(Gamma))
  Psi <- array(list(), dim=c(DIM, DIM))
  m <- apply(W, 1, sum)
  for(i in 1:DIM) {
    Gammai <- Gamma / m[i]
    Psi[i, i] <- list(solve(Gammai))
    if(i < DIM) for(j in (i+1):DIM) {
      if(W[i,j] == 1) {
        Lambdaij <- A %*% B %*% solve(A) / m[i]
        Psiij <- - solve(Gammai) %*% Lambdaij
        Psi[i,j] <- list(Psiij)
        Psi[j,i] <- list(t(Psiij))
      }
      else
        Psi[i,j] <- Psi[j,i] <- list(matrix(0, P, P))
    }
  }
  return(Psi)
}
Psi <- makePrecision(W, Gamma, B)
psi <- asMatrix(Psi)

phi <- chol(solve(psi)) %*% rnorm(P*DIM)
write.table(phi, file="phi.dat", row.names=FALSE, col.names=FALSE)

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

