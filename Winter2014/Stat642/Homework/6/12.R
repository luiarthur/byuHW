# Power is the probability of rejecting H0 when H1 is true
rm(list=ls())

N <- c(1,4,16,64,100)
a <- .05

plot.power <- function(my.beta,legend.place,fr,to,main){

  for (i in 1:length(N)) {
    add <- ifelse(i==1,F,T)
    curve(my.beta(x,N[i]),from=fr,to=to,col=i,add=add,lwd=2,
          xlab=expression(mu),ylab=expression(beta(mu)),
          main=main)
  }
  legend(legend.place,col=1:length(N),lwd=3,
         legend=paste("Power curve for n =",N))
  abline(h=.05)
}
#8.12a: ####################################################

my.beta.a <- function(m,n,m.0=0,s=1,a=.05){
  c <- qnorm(1-a)
  1 - pnorm(c+(m.0-m)/(s/sqrt(n)))
}

#plot.power(my.beta.a,"right",0,5)
###########################################################

#8.12b: ####################################################

my.beta.b <- function(m,n,m.0=0,s=1,a=.05){
  # Everything but the stuff in the middle
  c <- qnorm(1-a/2)
  1 - (pnorm(c+(m.0-m)/(s/sqrt(n))) - pnorm(-c+(m.0-m)/(s/sqrt(n))))
}

#plot.power(my.beta.b,"left",-10,5)
###########################################################

plot.12a <- function(){
  plot.power(my.beta.a,"right",0,5,main="12a: Power Curve")
}
plot.12b <- function(){
  plot.power(my.beta.b,"left",-10,5,main="12b: Power Curve")
}
