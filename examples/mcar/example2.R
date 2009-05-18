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

phi.true <- chol(solve(psi)) %*% rnorm(P*DIM)
beta.true <- 1:3 - 2
X <- matrix(, DIM*P, P)
for(i in 1:DIM)
  X[(i-1)*P + 1:P, 1:P] <- diag(P)
mu.true <- exp(X %*% beta.true + phi.true)
y <- rpois(DIM*P, mu.true)
write.table(y, file="y.dat", row.names=FALSE, col.names=FALSE)

summary(glm(y ~ X-1, family=poisson))

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
plot(phi[,1:3])
vphi <- var(phi)
