#1a
z <- x^(a^b)
#1b
z <- (x^a)^b
#1c
z <- 3*x^3 + 2*x^2 + 6*x + 1
#1d
z <- floor(x*100)-floor(x*10)*10
#1e
z <- z + 1

#2a
c(1:8,7:1)
#2b
rep(1:5,1:5)
#2c
matrix(1,3,3)-diag(c(1,1,1))
#2d
matrix(c(0,2,3,0,5,0,7,0,0),nrow=3,byrow=T)

#3
vec <- c(1,1)
polar <- function(vec){
  r <- sqrt( vec[1]^2 + vec[2]^2 )
  theta <- atan( abs(vec[2]/vec[1]) )
  if ( (vec[1] >= 0) & (vec[2] >= 0) ) {
    c( r, theta)
  } else if ( (vec[1] < 0) & (vec[2] >=  0) ) {   
    c(r, pi/2+theta) 
  } else if( (vec[1] < 0) & (vec[2] <  0) ) { 
    c(r, pi+theta)
  } else {
    c (r, 1.5*pi+theta) 
  }
}

polar(vec) # vec is the vector c(x,y)

#4
subset( 1:100,(1:100%%2!=0)&(1:100%%3!=0)&(1:100%%7!=0) )

#5
queue <- c("Steve", "Russell", "Alison", "Liam")
#5a
queue <- c(queue, "Barry")
#5b
queue <- queue[-1]
#5c
queue <- c("Pam", queue)
#5d
queue <- queue[queue != "Barry"]
#5e
queue <- queue[queue != "Alison"]

which(queue == "Russell")
