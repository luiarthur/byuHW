
x <- c('Sunday'=1,
       'Monday'=5,
       'Tuesday'=9,
       'Wednesday'=6,
       'Thursday'=5,
       'Friday'=6,
       'Saturday'=0)

#1
n <- length(x)
E <- rep(sum(x) / n,n)
Xsq <- sum((x-E)^2 / E)  
Xsq

#2
pchisq(Xsq,n-1,lower.tail=F)

#3

my.p.chiSq <- function(x,N=50000){

  n <- length(x)
  E <- rep(sum(x)/n, n)
  Xsq <- sum((x-E)^2 / E)

  one.it <- function(dummy) {
    samp <- rmultinom(1,sum(x), rep(1/n,n))
    n.samp <- length(samp)
    E.samp <- rep(sum(samp)/n.samp,n.samp)
    sum((samp-E.samp)^2 / E.samp)
  }
  
  Xsq.samp <- apply(matrix(0,N,1),1,one.it)
  
  mc.pval <- mean(Xsq.samp>=Xsq)
  CI <- qnorm(c(.025,.975),mc.pval,sqrt(mc.pval*(1-mc.pval)/N))
  list("MC.pval"=mc.pval,"CI"=CI)
}
my.p.chiSq(x)

#4: Put in layman's terms
#We reject the null hypothesis that local birthdays are 
#equally likely on all days of the week at the .05 confidence level.

#5
#The Monte Carlo approach can be taken when there is 
#no convenient or closed form solution to computing 
#a test statistic.
#However, the Monte Carlo approach takes simulation
#time and computing resources. So it may not be the 
#best choice when a closed form solution is known.
#Moreover, Monte Carlo simulations are approximations
#at best. Closed form solutions are accurate.
#DBD: Closed form solutions may be accurate, but not necessarily.
#An asymptotical closed-form solution *may be* very poor for
#a finite sample size.  Monte Carlo methods are not
#"approximations at best"; they are exact solutions up to
#Monte Carlo error, but you can control Monte Carlo error.
#On the other hand, it is difficult to control or quantity
#the amount of error when applying an asymptotic distribution
#when the sample size is finite.

#6
prob <- c(.05,.2,.2,.2,.2,.1,.05)
power.sim <- function(x,prob,N=10000){
  n <- length(x)

  Xsq <- NULL
  Xsq.samp <- NULL

  for (i in 1:N){
    obs <- rmultinom(1,sum(x),rep(1,n))
    Xsq[i] <-sum((obs-E)^2 / E)
    samp <- rmultinom(1,sum(x),prob)
    Xsq.samp[i] <- sum((samp-E)^2 / E)
  }
  
  crit <- quantile(Xsq,.95)

  den1 <- density(Xsq.samp)
  den2 <- density(Xsq)

  pdf('plot.pdf')
    plot(0,0,cex=0,ylim=c(0,.14),xlim=c(0,40),
      xlab=expression(chi^2),ylab='Density',
      main=expression('Density of the Sampling Distribution of the Test Statistic'))

    polygon(c(0, den1$x[den1$x>0 & den1$x < crit], crit), 
            c(0, den1$y[den1$x>=0 & den1$x <= crit], 0),col="red")
    polygon(c(crit, den2$x[den2$x>crit & den2$x < 40], 40), 
            c(0, den2$y[den2$x>=crit & den2$x <= 40], 0),col="blue")

    lines(den1,col='red',lwd=3)
    lines(den2,col='blue',lwd=3)
    curve(dchisq(x,n-1),from=0,to=50,lwd=3,col='gold',add=T)

    legend(20,.12,legend=c('Ho: Asymptotic','Ho: Finite Sampling','Ha'),
           col=c('gold','blue','red'),lwd=3)

    text(7,.02,'Type II Error',font=2)
    text(17,.01,'Type I Error',font=2)
  dev.off()

  mc.power <- mean(Xsq.samp >= quantile(Xsq,.95))
  CI.mc <- mc.power + c(-1,1) * qnorm(.975) * sqrt(mc.power*(1-mc.power)/N)

  asym.power <- mean(Xsq.samp >= qchisq(.95,n-1))
  CI <- asym.power + c(-1,1) * qnorm(.975) * sqrt(asym.power*(1-asym.power)/N)

  list(c("MC.power"=mc.power, "CI.lower"=CI.mc[1], "CI.upper"=CI.mc[2]), 
       c("Asymptotic.Power"=asym.power, "CI.lower"=CI[1], "CI.upper"=CI[2]))
}
power.sim(x,prob)
