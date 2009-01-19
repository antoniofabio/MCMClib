Fn <- function(X0, Y) {
  writeLines(as.character(X0), f <- file("tmp_X0.dat"))
  close(f)
  as.numeric(system("./distrCompute 1 tmp_X0.dat", input=as.character(Y), intern=TRUE))
}

X0 <- c(-1.0, 0.0, 1.0)
Y <- c(-0.5, 0.5, 0.0)
Fn(X0, Y)
#[1] 0.000000 0.666667 1.000000

X0 <- rnorm(5000)
Y <- rnorm(100000)
system.time(F <- Fn(X0, Y))
hist(F, breaks=30)

FnDist <- function(X0, F, X) {
  M <- dim(X)[3]
  writeLines(as.character(F), f <- file("tmp_F.dat")); close(f)
  writeLines(as.character(X0), f <- file("tmp_X0.dat")); close(f)
  as.numeric(system(paste("./distrDist", M, "tmp_F.dat tmp_X0.dat"),
         intern=TRUE, input=as.character(X)))
}
X <- rnorm(10000)
system.time(D <- FnDist(X0, F, array(X, dim=c(10000, 1, 1))))
plot(D, type="l", log="y")
