#i)

p <- function(xbar,s=1,n=9) {
  2 * (1 - pnorm(abs(xbar)/(s/sqrt(n))))
} 

po <- function(x,n=9){

  m <- function(x){
    .5 * sqrt(n/(2*pi)) * exp(-n*x^2/2) +
    .5 * 1/sqrt(2*pi*(1/n+1)) * exp(-x/(2*(1/n+1)))
  }
 
  f <- function(x){
    sqrt(n/2/pi) * exp(-n*x^2/2)
  }

  f(x)/m(x)/2
}

po.n <- function(n) po(qnorm(.05/2)/sqrt(n),n)

plot.d1 <- function(){
  curve(p(x),from=0,to=1.5,xlab=expression(bar(x)),col="red",lwd=3) # pval
  curve(po(x),from=0,to=1.5,xlab=expression(bar(x)),col="blue",lwd=3,add=T) # post
  title(expression(paste("Probability vs. ",bar(x))))
  legend("topright",legend=c("P-value","Bayes"),col=c("red","blue"),lwd=3)
  abline(v=.104)
  text(.25,.9,expression(bar(x) %~~% .104))
}

plot.d2 <- function(N=10^7){
  curve(po.n(n),ylim=c(0,1),from=10^-10,to=N,col="blue",lwd=3,xname="n") #bayes
  lines(c(0,N),c(.05,.05),type='l',col='red',lwd=3) #pval
  title(expression(paste("Probability vs. ",bar(x))))
  legend("right",legend=c("P-value","Bayes"),col=c("red","blue"),lwd=3)
  abline(h=1)
}

pdf("d1.pdf")
  plot.d1()
dev.off()  
pdf("d2.pdf")
  plot.d2()
dev.off()  
