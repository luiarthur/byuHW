#22a:
c <- 0:10

f <- function(c,n=10,p=.5){
  choose(n,c) * p^c * (1-p)^(n-c)
}

f(c)

#c:
sum(f(0:2))
# => c = 2

#:power
f(c,10,.25)
sum(f(0:2,10,.25))

#22b:
f(c,10,.5)
(size <- sum(f(6:10,10,.5)))

power <- function(p) {
  sum <- 0
  for (i in 6:10) {
    sum <- sum + f(i,10,p)
  }
  sum
}  

plot.22 <- function(){
  curve(power(x),from=0,to=1,col='red',lwd=3,
        xlab=expression(p),
        ylab=expression(beta(p)),
        main="22: Power Curve")
  curve(power(x),from=0,to=.5,col='blue',lwd=3,add=T)
  lines(c(0,.5,.5),c(size,size,0))
  text(.2,.4,paste("size =",size))
  legend("topleft",lwd=3,col="red",legend="Power Curve")
}
