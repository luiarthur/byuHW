# y_i ~ POI(exp(beta*x_i))
#rm(list=ls())
options(width=as.integer(120))

beta <- seq(0,.8,by=.05)

x <- c(.91,2.34,2.54,1.57,1.63,1.08,.3,.9,.08,2.01,1.5,1.2)

generate.y <- function(x,b){
  y <- NULL
  for (i in 1:length(x)){
    y[i] <- rpois(1,exp(b*x[i]))
  }  
  y
}

#Bayes #############################################################
est.mode <- function(x){
  d <- density(x)
  d$x[which.max(d$y)]
}
bayesEstimate <- function(N=11000,burn=1000,cs,prior=c('normal','uniform')[1]) {
  #B <- matrix(0,length(beta),5)
  #colnames(B) <- c('Acceptance','Mean.of.Draws','True.Beta','Bias','MSE')
  B <- NULL

  for (z in 1:length(beta)){
    b <- beta[z]
    y <- generate.y(x,b)

    out <- rnorm(1,.4,.5)
    log.target <- function(b) sum(-exp(x*b) + x*y*b) - ((b-.4)/.5)^2/2
    if (prior=='uniform'){
      out <- runif(1,-2,.2)
      log.target <- function(b) sum(-exp(x*b) + x*y*b)
    }

    cnt <- 0

    for (i in 2:N){
      out[i] <- out[i-1]
      cand <- rnorm(1,out[i],cs[z])
      if (!((prior=='uniform') & (abs(cand)>2))) {
        if ( (log.target(cand) - log.target(out[i])) > log(runif(1)) ){
          out[i] <- cand
          cnt <- cnt + 1
        }
      }
    }

    out <- out[-(1:burn)]
    #plot(out,type='l')
    accpt <- cnt/N
    #B[z,] <- matrix(c(accpt,mean(out),b,
    #                  mean(out)-b,var(out)+(mean(out)-b)^2),nrow=1)
    #B[z] <- ifelse(prior=='normal',mean(out),out[which.max(out)])
    B[z] <- ifelse(prior=='normal',mean(out),est.mode(out))

  }
  
  # output is length 17 vector
  B #plot(B[,4],pch=20); abline(h=0)
}

library(foreach)
library(doMC)
registerDoMC(16)

engine.normal <- function(){
  normalBayEst <- bayesEstimate(cs=c(rep(.5,2),rep(.7,8),rep(.2,7)),prior='normal')
  matrix(normalBayEst,nrow=1)
}
engine.uniform <- function(){
  unifBayEst <- bayesEstimate(cs=c(rep(.5,2),rep(.7,8),rep(.2,7)),prior='uniform')
  matrix(unifBayEst,nrow=1)
}

N <- 10000
bayes.Ests.norm <- foreach(k=1:N,.combine=rbind) %dopar% engine.normal()
bayes.Ests.unif <- foreach(k=1:N,.combine=rbind) %dopar% engine.uniform()

bias.norm <- apply(bayes.Ests.norm,2,mean) - beta 
MSE.norm <-apply(bayes.Ests.norm,2,var) + bias.norm^2 
results.norm <- cbind(bias.norm,MSE.norm)
colnames(results.norm) <- c('Bias','MSE')
#results.norm

bias.unif <- apply(bayes.Ests.unif,2,mean) - beta
MSE.unif <-apply(bayes.Ests.unif,2,var) + bias.unif^2
results.unif <- cbind(bias.unif,MSE.unif)
colnames(results.unif) <- c('Bias','MSE')
#results.unif

#MLE###############################################################
calc.mle.bias <- function(beta,R=1000){
  mle.results <- NULL
  for (i in 1:R){
    b <- beta
    y <- generate.y(x,b)
    f1 <- function(beta) sum(x * (-exp(x*beta) + y))
    f2 <- function(beta) sum(-x^2 * exp(x*beta))

    newton.raphson <- function(f1,f2,init=0,eps=.0000001){
      b <- init 
      while( abs(f1(b)) > eps ){
        b <- b - f1(b)/f2(b)
      }
      b
    }

    mle.results[i] <- newton.raphson(f1,f2,1)
  }
  bias <- mean(mle.results) - b
  mse  <- var(mle.results) + bias^2
  out <- matrix(c(bias,mse),nrow=1)
  colnames(out) <- c('Bias','MSE')
  out
}

###############For Non-Parallel ######################
#results.mle <- matrix(0,length(beta),2)             #
#colnames(results.mle) <- c('Bias','MSE')            #
#for (i in 1:length(beta)){                          #
#  results.mle[i,] <- calc.mle.bias(beta[i],R=10000) #
#}####################################################

#Parallel
results.mle <- foreach(k=1:length(beta),.combine=rbind) %dopar% 
                 calc.mle.bias(beta[k],10000)
#results.mle

results <- cbind(beta,results.mle,results.norm,results.unif)
colnames(results) <- c('True.Beta','MLE.Bias','MLE.MSE',      # red
                       'Bayes.Norm.Bias','Bayes.Norm.MSE',    # blue
                       'Bayes.Unif.Bias','Bayes.Unif.MSE')    # green
print(results)

# Plots
pdf("bias.pdf")
  ml <- max(results[,c(2,4,6)]) + .01
  ylim <- c(min(results[,c(2,4,6)]),ml)
  plot(results[,1],results[,2],ylim=ylim,col='red',
       main=expression(paste('Bias vs. ',beta)),
       xlab=expression(beta),ylab='Bias',pch=20,type='o')
  points(results[,1],results[,4],ylim=ylim,col='blue',pch=20,type='o')
  points(results[,1],results[,6],ylim=ylim,col='green',pch=20,type='o')
  legend("topright",legend=c('MLE','Bayes Est. Normal Beta',
         'Bayes Est. Uniform Beta'), col=c('red','blue','green'),lwd=3)
dev.off()

pdf("mse.pdf")
  ml <- max(results[,c(3,5,7)]) 
  ylim <- c(min(results[,c(3,5,7)]),ml)
  plot(results[,1],results[,3],ylim=ylim,col='red',
       main=expression(paste('MSE vs. ',beta)),
       xlab=expression(beta),ylab='MSE',pch=20,type='o')
  points(results[,1],results[,5],ylim=ylim,col='blue',pch=20,type='o')
  points(results[,1],results[,7],ylim=ylim,col='green',pch=20,type='o')
  legend("topright",legend=c('MLE','Bayes Est. Normal Beta',
         'Bayes Est. Uniform Beta'), col=c('red','blue','green'),lwd=3)
dev.off()

