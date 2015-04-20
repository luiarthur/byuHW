source("rfunctions.R")
n <- 10
A <- matrix(c(3,7,11))
Z <- diag(1,3) %x% matrix(rep(1,n))
E <- rnorm(n*3)
Y <- Z%*%A+E
