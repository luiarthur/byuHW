c <- 1.68377

b2 <- function(t) {
  if (t <= c/2 - 1){
    0
  } else if ( t <= (c-1)/2 ) {
    (2*t+2-c)^2/2
  } else if( t <= c/2 ) {
    1 - (c-2*t)^2/2
  } else {
    1
  }
}

b1 <- function(t){
  if (t <= -.05){
    0
  } else if (t <= .95) {
    t + .05
  } else {
    1
  }  
}

t <- seq(-.05,1,length=100)
b1t <- NULL
b2t <- NULL

for (i in 1:length(t)){
  b1t[i] <- b1(t[i])
  b2t[i] <- b2(t[i])
}

plot.13 <- function(){  
  plot(t,b2t,type='l',col='red',lwd=3,
       ylab=expression(beta(theta)),
       xlab=expression(theta),
       main="8.13b: Power Curves")

  lines(t,b1t,type='l',col='blue',lwd=3)
   
  legend("bottomright",legend=paste("Power of Test",1:2),col=c("blue","red"),lwd=3)
}  
