options(digits=15)

#2
h2 <- function(x=6.6,n=8){
  sum <- 0
  for (i in 0:n){
    sum <- sum + x^i
  }
  sum
}

#3
h3 <- function(x=6.6,n=8){
  (1- x^(n+1)) / (1-x)
}

#4
h4a <- function(x=6.6,n=8){
  sum <- 0
  i <- 0
  while ( i <= n ){
    sum <- sum + x ^ i 
    i <- i + 1
  }
  sum
}

h4b <- function(x=6.6,n=8){
  sum(rep(x,n+1)^(0:n))
}

#6
h6a <- function(x= (1:500) ){
  n <- length(x)
  logProd <- 0
  for (i in 1:n){
    logProd <- logProd + log(x[i])
  }
  exp( 1/n * logProd)
}

h6b <- function(x= (1:500) ){
  n <- length(x)
  exp (1/n * sum(log(x)))
}

#7
q7 <- function(v=(2:30 * 3)){
  sum( v[ 1:floor(length(v)/3) * 3 ] )
}

# FUNCTIONS FOR JARED:

h2()     #input: (x,n)
h3()     #input: (x,n)
h4a()    #input: (x,n)
h4b()    #input: (x,n)
h6a()    #input: a vector
h6b()    #input: a vector
q7()     #input: a vector

