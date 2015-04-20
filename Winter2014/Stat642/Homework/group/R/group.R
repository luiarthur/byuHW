dat <- read.csv("../Data/jazz.csv")

datT <- dat[dat$strike==T,]
datF <- dat[dat$strike==F,]


hist(datF[,1],freq=F)
curve(dgamma(x,6,scale=6),from=0,60,add=T)

up.Lam <- function(x,a,b){
  n <- length(x)
  new.a <- a+sum(x)
  new.b <- 1/(1/b+n)
  list("a"=new.a,"b"=new.b)
}

#gibbs <- function(x=dat,init=6,B=1000){
#  lam <- NULL
#  lam[1] <- init
#  a <- 6; b <- 6
#  UP <- matrix(0,B,2)
#  UP[1,1] <- a; UP[1,2] <- b
#  for (i in 2:B) {
#    up <- up.Lam(x,a,b)
#    UP[i,1] <- up$a; UP[i,2] <- up$b
#    lam[i] <- rgamma(1,up$a,rate=up$b)
#  }
#  list(lam,UP)
#}
#
#outT <- gibbs(x=datT[,1],B=10000)[[1]]
#outF <- gibbs(x=datF[,1],B=10000)[[1]]

lamT <- up.Lam(datT[,1],6,6)
lamF <- up.Lam(datF[,1],6,6)
outT <- rgamma(100000,lamT$a,scale=lamT$b)
outF <- rgamma(100000,lamF$a,scale=lamF$b)

diff <- outT - outF
plot(density(outT),xlim=c(30,40))
lines(density(outF))

dis <- NULL
B <- 10000
perc <- seq(0,.05,length=B)
for (i in 1:B){
  l <- perc[i]
  u <- l + .95
  dis[i] <- dist(quantile(diff,c(l,u)))
}  

ind <- which.min(dis)
lower.perc <- perc[ind]
upper.perc <- lower.perc + .95
dis[ind]

#1: HPD
plot.hpd <- function() {
  plot(density(diff),main=expression(paste("95% HPD for the Difference in ",lambda)),lwd=3,col="red")
  hpd <- c(round(quantile(diff,c(lower.perc,upper.perc)),4))
  abline(v=hpd)
  abline(h=.058)
  legend("topright",legend=c(paste("lower =",hpd[1],4), paste("upper =",hpd[2],4)))
  text(-4,.07,".058")
}

#2: Equal-tailed Credible Interval
plot.eq <- function() {
  plot(density(diff),main=expression(paste("95% Equal-Tailed Credible Interval for the Difference in ",lambda)),lwd=3,col="blue")
  eqCred <- c(round(quantile(diff,c(.025,.975)),4))
  abline(v=eqCred)
  legend("topright",legend=c(paste("lower =",eqCred[1]),paste("upper =",eqCred[2])))
}

#3: Confidence Interval
plot.ci <- function() {
  plot(density(diff),main=expression(paste("95% Confidence Interval for the Difference in ",lambda)),lwd=3,col="orange")
  CI <- round(qnorm(c(.025,.975),mean(diff),sd(diff)),4)
  names(CI) <- c("2.5%","97.5%")
  abline(v=CI)
  legend("topright",legend=c(paste("lower =",CI[1]),paste("upper =",CI[2])))
}

pdf("hpd.pdf"); plot.hpd(); dev.off()
pdf("eq.pdf");  plot.eq(); dev.off()
pdf("ci.pdf");  plot.ci(); dev.off()


#1: HPD
hpd
#2: Equal-tailed Credible Interval
eqCred
#3: Confidence Interval
CI

