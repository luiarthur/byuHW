rm(list=ls())
soil <- read.csv("soil.csv")[,-1]
#plot(soil,pch=20,col='brown',main="SWC Vs. CWSI")

rmvn <- function(n=1,mu=0,Sigma=1){
  draws <- mu + t(chol(Sigma)) %*% rnorm(n)
  draws
}

#install.packages("geoR")
#install.packages("LatticeKrig")
library(geoR)
library(LatticeKrig)      #Load rdist function
# small nu  => crooked, jagged line
# small phi => small seasonal effects

                                    # var(Y)       #range(X)   
#GP <- function(data=soil,nu=2,K=101,s2.start.val=1,phi.start.val=1,pred=data[1,],plot=F){
GP <- function(data=soil,nu=2,K=101,s2.start.val=.5*var(soil[,2]),phi.start.val=.001,pred=data[1,],plot=F,export.plot=F){
  N <- nrow(data)
  SWC <- data$SWC
  CWSI <- data$CWSI

  obs <- as.geodata(cbind(SWC,CWSI,rep(0,N)),data.col=1,coords.col=2:3)
  gp.fit <- likfit(obs,cov.model="matern",kappa=nu,fix.kappa=TRUE,
                   ini.cov.pars=c(s2.start.val,phi.start.val), trend="cte")
  phi <- 1/gp.fit$phi
  s2 <- gp.fit$sigmasq
  mu <- gp.fit$beta
  tau2 <- gp.fit$tausq

  #pred.seq <- seq(min(CWSI),max(CWSI),length=K)
  #pred.seq <- c(seq(min(CWSI),max(CWSI),length=K-1),pred[1])
  pred.seq <- c(seq(0,1,length=K-1),pred[1])
  D <- rdist(c(pred.seq,CWSI))
  V <- s2*Matern(D,alpha=phi,nu=nu) ##V = Sigma_Y
  EV <- mu + V[1:K,K+(1:N)] %*% solve(V[K+(1:N),K+(1:N)]+tau2*diag(N)) %*% (SWC-mu)
  cond.Var <- diag((V[1:K,1:K]+tau2*diag(K))-V[1:K,K+(1:N)] %*% 
              solve(V[K+(1:N),K+(1:N)] + tau2*diag(N))%*%t(V[1:K,K+(1:N)]))

  upper <- qnorm(0.975,mean=EV,sd=sqrt(cond.Var))
  lower <- qnorm(0.025,mean=EV,sd=sqrt(cond.Var))

  if (plot) {
    if (export.plot) pdf("predict.pdf")
    plot(pred.seq[-K],EV[-K],type="l",lwd=3,xlab="CWSI",ylab="SWC", #####
         #xlim=range(CWSI), ylim=c(20,29),col="red",
         xlim=range(pred.seq), ylim=c(min(lower),max(upper)),col="red",
         main=paste("Gaussian Process Predictions, ",
                    expression(nu), "=",round(nu,3)))
    points(CWSI[-K],SWC[-K],pch=19,cex=0.5)
    lines(pred.seq[-K],lower[-K],col="blue")
    lines(pred.seq[-K],upper[-K],col="blue")
    legend("topright",legend=c("Prediction Estimate","95% Prediction Bands"),
           col= c("red","blue"),lwd=2)
    if (export.plot) dev.off()
  }  
  #D <- rdist(CWSI)
  #V <- s2*Matern(D,alpha=phi,nu=nu) + tau2*diag(N)
  #X <- cbind(rep(1,N))
  #Y <- cbind(SWC)
  #b.var <- solve(t(X) %*% solve(V) %*% X)
  #b <- b.var %*% t(X) %*% solve(V) %*% Y
 
  pred.in <- ifelse(lower[K] < pred[2] & pred[2] < upper[K],T,F)

  list("pred.in"=pred.in,"gp.fit"=gp.fit,"pred"=EV[K],
       "up"=upper[-K],"lo"=lower[-K],"ev"=EV[-K])
}

CV <- function(reg=2,data=soil){

  library(doMC)
  library(foreach)
  registerDoMC(reg)

  n <- nrow(data)

  leave.1.out <- function(i){
    estimates <- GP(data[-i,],pred=data[i,],plot=F)
    estimates$pred.in
  }

  coverage <- foreach(i=1:n,.combine=cbind) %dopar% leave.1.out(i)
  coverage
}

cv <- CV(16)
n <- length(cv)
p <- mean(cv)
#coverage.CI <- p + c(-1,1)*qnorm(.975)*sqrt((p*(1-p)/n)) # same as below
coverage.CI <- qnorm(c(.025,.975),p,sqrt((p*(1-p)/n)))
coverage.CI[2] <- min(1,coverage.CI[2])
coverage <- cbind(p,t(coverage.CI))
colnames(coverage) <- c("Est.Coverage","CI.lo","CI.hi")

est <- GP(nu=2,plot=T,export.plot=T)
beta <- est$gp.fit$beta
se <- sqrt(est$gp.fit$beta.var)
beta.CI <- qnorm(c(.025,.975),beta,se)

plot.resid <- function(){
  par(mfrow=c(2,1))
  resids <- est$gp.fit$model.components$residuals
  #plot(resids,pch=20); abline(h=0)
  qqnorm(resids,pch=20)
  hist(resids,col='green')
  par(mfrow=c(1,1))
}

widths <- est$up - est$lo
min.widths <- min(widths)
max.widths <- max(widths)
which.min.widths <- which.min(widths) / 100
which.max.widths <- which.max(widths) / 100
min.ev <- est$ev[which.min(widths)]
max.ev <- est$ev[which.max(widths)]

set <- seq(1,100,length=20) 
pred.table <- cbind(set/100,est$ev[set],est$lo[set],est$up[set])
colnames(pred.table) <- c("CWSI","Est.SWC","PI.Lo","PI.Up")
