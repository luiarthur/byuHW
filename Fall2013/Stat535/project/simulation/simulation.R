
# MODEL MISSPECIFICATION: Take 2
#
# True Model:
# Y = b0+b1X1+b2(X1)^2 + eps
#
# Misspecified Model:
# Y = b0 + b1X1 + eps


simulate <- function(b,n=100,N=100000){
  #b <- b[1,]
  b0 <- b[1]; b1 <- b[2]; b2 <- b[3]

  lm.F <- function(y){
    cf <- matrix(c(lm(y ~ x)$coefficients,0),ncol=1)
    rownames(cf) <-c('b0','b1','b2')
    cf
  }
  lm.T <- function(y) {
    cf <- matrix(lm(y ~ x + I(x^2))$coefficient,ncol=1)
    rownames(cf) <- c('b0','b1','b2')
    cf
  }

  # N = number of iterations for simulation
  # n = number of data points
  eps <- matrix(rnorm(n*N),n,N)
  x <- rnorm(n,1,1)
  Y <- b0 + b1*x + b2*x^2 + eps

  mod.F <- apply(Y,2,lm.F) #lm(Y[,1]~x)
  #mod.T <- apply(Y,2,lm.T) #lm(Y[,1]~x+I(x^2))
  #plot(x,Y[,1]); abline(mod.F)
  #bh.T <- mod.T$coefficients
  #xx <- seq(-2,4,length=1000)
  #yy <- bh.T[1] + bh.T[2]*xx + bh.T[3]*xx^2
  #lines(xx,yy)
  my.lm <- function(y){
    CI <- confint(lm(y~x))
    b0.in <- ifelse((CI[1,1] < b0) & (b0 < CI[1,2]),T,F)
    b1.in <- ifelse((CI[2,1] < b1) & (b1 < CI[2,2]),T,F)
    b2.in <- ifelse(b2==0,T,F)
    matrix(c(b0.in,b1.in,b2.in),ncol=1)
  }

  mod <- apply(Y,2,my.lm)
  b0.in <- mean(mod[1,])
  b1.in <- mean(mod[2,])
  b2.in <- mean(mod[3,])
  in.CI <- matrix(c(b0.in,b1.in,b2.in),nrow=1)

  bias.b0 <- mean(mod.F[1,]) - b0
  bias.b1 <- mean(mod.F[2,]) - b1
  bias.b2 <- mean(mod.F[3,]) - b2

  mse.b0 <- var(mod.F[1,]) + bias.b0^2
  mse.b1 <- var(mod.F[2,]) + bias.b1^2
  mse.b2 <- var(mod.F[3,]) + bias.b2^2

  bias <- cbind(bias.b0, bias.b1, bias.b2)
  mse <- cbind(mse.b0, mse.b1, mse.b2)

  true <- c(b0,b1,b2)
  results <- rbind(true,bias,mse,in.CI)
  rownames(results) <- c('True Stats','Bias','MSE','in.CI%')
  colnames(results) <- c('b0.hat','b1.hat','b2.hat')
  results
}

library(foreach)
library(doMC)
registerDoMC(16)
useParallel <- TRUE

b0 <- seq(0,.5,length=100)#c(1,1,4)
b1 <- seq(0,.5,length=100)#c(1,1,4)
b2 <- seq(0,.5,length=100)#c(0,.05,1)
b <- cbind(b0,b1,b2)

results <- list()
if (useParallel){
  sim.time <- system.time(
    results <- foreach(k=1:nrow(b)) %dopar% simulate(b[k,],N=1000)
  )
  print(sim.time)
} else {
  begin.time <- Sys.time()
  results <- list() 
  for (k in 1:nrow(b)){
    results[[k]] <- simulate(b[k,])
  }
  print(Sys.time() - begin.time)
}

bias.b0 <- NULL
bias.b1 <- NULL
bias.b2 <- NULL
MSE.b0 <- NULL
MSE.b1 <- NULL
MSE.b2 <- NULL
cvg.b0 <- NULL
cvg.b1 <- NULL
cvg.b2 <- NULL
for(i in 1:length(results)){
  bias.b0[i] <- results[[i]][2,1]
  bias.b1[i] <- results[[i]][2,2]
  bias.b2[i] <- results[[i]][2,3]
  MSE.b0[i] <- results[[i]][3,1]
  MSE.b1[i] <- results[[i]][3,2]
  MSE.b2[i] <- results[[i]][3,3]
  cvg.b0[i] <- results[[i]][4,1]
  cvg.b1[i] <- results[[i]][4,2]
  cvg.b2[i] <- results[[i]][4,3]
}

par(mfrow=c(3,1))

plot(b0,bias.b0,ty='o',pch=20,col='red',
     main=expression(paste('BIAS Vs. ',beta)),
     xlab=expression(paste('(',beta[1],beta[2],beta[3],')')),
     ylim=c(min(bias.b0,bias.b1,bias.b2),max(bias.b0,bias.b1,bias.b2)))
lines(b1,bias.b1,ty='o',pch=20,col='blue',
      main=expression(paste('BIAS Vs. ',beta[2])))
lines(b2,bias.b2,ty='o',pch=20,col='green',
      main=expression(paste('BIAS Vs. ',beta[3])))

plot(b0,MSE.b0,ty='o',pch=20,col='red',
     main=expression(paste('MSE Vs. ',beta)),
     xlab=expression(paste('(',beta[1],beta[2],beta[3],')')),
     ylim=c(min(MSE.b0,MSE.b1,MSE.b2),max(MSE.b0,MSE.b1,MSE.b2)))
lines(b1,MSE.b1,ty='o',pch=20,col='blue',main=expression(paste('MSE Vs. ',beta[2])))
lines(b2,MSE.b2,ty='o',pch=20,col='green',main=expression(paste('MSE Vs. ',beta[3])))

plot(b0,cvg.b0,ty='o',pch=20,col='red',
     main=expression(paste('Coverage Vs. ',beta)),
     xlab=expression(paste('(',beta[1],beta[2],beta[3],')')),
     ylim=c(min(cvg.b0,cvg.b1,cvg.b2),max(cvg.b0,cvg.b1,cvg.b2)))
lines(b1,cvg.b1,ty='o',pch=20,col='blue',
      main=expression(paste('Coverage Vs. ',beta[2])))
lines(b2,cvg.b2,ty='o',pch=20,col='green',
      main=expression(paste('Coverage Vs. ',beta[3])))

results[[1]]
results[[10]]
results[[20]]
