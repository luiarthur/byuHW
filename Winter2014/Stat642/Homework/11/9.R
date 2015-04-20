x <- c(10,7,8,13,8,9,5,7,6,8,3,6,6,3,5)
y <- sum(x)

P <- function(x,t) {
  y <- sum(x)
  n <- length(x)
  choose(y,t) * (n-1)^(y-t) / n^y
}
