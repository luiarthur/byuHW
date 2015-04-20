x <- 1:20
y <- c(1,6,16,23,27,39,31,30,43,51,63,70,88,97,91,104,110,113,149,159)

# Plots: ###################################################################
par(mfrow=c(1,2),cex.main=.9)

# 1a)
plot1a <- function() {
  plot(x,y,pch=20,cex=1.5,main="Number of Cases Vs. Time Period",
       xlab="Time Period",ylab="Number of Cases",cex.main=.9)
}
plot1a()

# 1b)
plot1b <- function() {
  plot(log(x),log(y),pch=20,cex=1.5,cex.main=.9,main="Log Number of Cases Vs. 
       \n Log Time Period", xlab="Log Time Period",ylab="Log Number of Cases")
}
plot1b()

plotall <- function() {
  par(mfrow=c(1,2))
    plot1a()
    plot1b()
  par(mfrow=c(1,1))
}
plotall()
############################################################################

# 1c)
iwls <- function(X,y,b=NULL,eps=1e-5,maxits=500) {
  X <- cbind(1,X)
  if (is.null(b)) b <- matrix(0,ncol(X),1)
  diff <- 1
  its <- 0

  while(diff>eps && its<maxits) {
    old.b <- b
    Xb <- X%*%b
    W <- diag(as.numeric(exp(Xb)))
    z <- Xb + y/exp(Xb) - 1
    tW <- t(X)%*%W
    b <- solve(tW%*%X,tW %*% z)
    diff <- sqrt(sum((b-old.b)^2))
    its <- its + 1
  }

  rownames(b) <- c("Intercept","Slope")
  #cat("Iterations: ",its,"\n")
  cov.b <- solve(tW%*%X)
  colnames(cov.b) <- rownames(cov.b) <- rownames(b)

  list("b"=b,"cov"=cov.b,"its"=its)
}

out <- iwls(log(x),y)
out$b

# 1d)
poismod <- glm(y~log(x),family=poisson)
summary(poismod)

# 1e)
out$cov
M <- matrix(0,2,3)
M[1,] <- c(out$b[1],qnorm(c(.025,.975),out$b[1],sqrt(diag(out$cov)[1])))
M[2,] <- c(out$b[2],qnorm(c(.025,.975),out$b[2],sqrt(diag(out$cov)[2])))
colnames(M) <- c("Estimate","Lower 95% CI","Upper 95% CI")
rownames(M) <- c("Intercept","Slope")
M

# Since neither confidence intervals contain 0, both betas are significant at
# the 95% confidence level.
