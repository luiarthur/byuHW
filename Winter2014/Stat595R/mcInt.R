f <- function(x,t) {
  1/t * exp(-x/t)
}

F <- function(x,t) {
  1- exp(-x/t)
}

S <- function(x,t) {
  exp(-x/t)
} 


N <- 10000000
ub <- 10000
u <- runif(N,0,ub)
mean(u*f(u,1) / (1/ub))
mean(S(u,1) / (1/ub))
mean(F(u,1) / (1/ub))

# So: 
  E[x] = Int_0^Inf{S(x)}
