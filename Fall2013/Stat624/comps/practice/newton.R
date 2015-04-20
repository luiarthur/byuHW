f <- function(x) x^2 - 2*x + 1
f1 <- function(x) 2*x - 2
f2 <- function(x) 2

eps <- 10^(-10)
x1 <- 100
x0 <- 1000
i <- 0

while(abs(x1-x0) > eps) {
  i <- i + 1
  print(i)
  x0 <- x1
  x1 <- x0 - f1(x0)/f2(x0)
  
}

x1
