display <- function(x,b) {
  CI <- matrix(c(mean(x),qnorm(c(.025,.975),mean(x),sd(x)/sqrt(b))),1)
  colnames(CI) <- c("Estimate","CI.Lower","CI.Upper")
  CI
}

a <- 1
b <- 1
f <- function(x) gamma(a+b) / (gamma(a)*gamma(b)) * x^(1-a) * (1-x)^(b-1)
h <- function(x) x*f(x)

mcInt <- function(h,c,d,B=10^4) {
  X <- runif(B,c,d)
  g <- 1/(d-c)
  draws <- h(X) / g
  display(draws,B)
}

mcInt(h,0,1,10^6)
