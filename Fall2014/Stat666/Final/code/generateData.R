a1 <- matrix(c(0,1,0,0,0,0,
               1,1,1,0,0,0,
               0,1,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0),6,6,byrow=T)

a2 <- matrix(c(0,0,0,1,1,1,
               0,0,0,1,0,1,
               0,0,0,1,1,1,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0),6,6,byrow=T)

a3 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               1,0,0,0,0,0,
               1,1,0,0,0,0,
               1,1,1,0,0,0),6,6,byrow=T)

a4 <- matrix(c(0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,0,0,0,
               0,0,0,1,1,1,
               0,0,0,0,1,0,
               0,0,0,0,1,0),6,6,byrow=T)

N <- 100
K <- 4
Z <- matrix(sample(0:1,N*K,replace=T),N,K)
constrain.z <- function(z) { # constrain each row of Z to have at least one class
  k <- length(z)
  while (sum(z)==0) {
    z <- sample(0:1,k,replace=T)
  }
  z
}
Z <- t(apply(Z,1,constrain.z))

A <- matrix(c(a1,a2,a3,a4),K,byrow=T)
ZA <- Z%*%A
sigX <- .5
E <- matrix(rnorm(prod(dim(ZA)),0,sigX),dim(ZA)[1],dim(ZA)[2])
Y <- ZA + E



