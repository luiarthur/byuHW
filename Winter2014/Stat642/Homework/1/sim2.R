# 7.2
x <- c(22.0, 23.9, 20.9, 23.8, 25.0, 24.0, 21.7, 23.8, 22.8, 23.1, 23.1, 23.5, 23.0, 23.0)

xbar <- mean(x)
sumlogx <- sum(log(x))
sumx <- sum(x)
n <- length(x)

f <- function(a) {
  -n*lgamma(a) - n*a*log(xbar/a) +(a-1)*sumlogx - sumx/(xbar/a)
}

# Maximize the log likelihood
optimize(f, c(0,1000), tol = 0.0001, maximum=T)


f1 <- function(a) {
  -n*digamma(a) - n*log(xbar/a) + sumlogx
}

# Solve for where the first derivative = 0
a <- uniroot(f1,c(.00001,1000))$root
#Two methods, same answer: a = 514.3354 => b = xbar/a = .0449401

# 7.10
x <- c(22.0, 23.9, 20.9, 23.8, 25.0, 24.0, 21.7, 23.8, 22.8, 23.1, 23.1, 23.5, 23.0, 23.0)

sumlogx <- sum(log(x))
n <- length(x)

b <- max(x)
a <- n / (n*log(b) - sumlogx)

