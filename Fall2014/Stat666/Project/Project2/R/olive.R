o2 <- as.matrix(read.table("oliver2a.dat",header=T)) 
o4 <- as.matrix(read.table("oliver4a.dat",header=T))
  
#Book data, for testing my EM algorithm. See p. 84 of book.
dat <- matrix(c(35,NA,2.80,
                35,4.9,NA,
                40,30,4.38,
                10,2.8,3.21,
                6,2.7,2.73,
                20,2.8,2.81,
                35,4.6,2.88,
                35,10.9,2.9,
                35,8,3.28,
                30,1.6,3.2),10,3,byrow=T)

EM <- function(X,thresh=1e-5,maxit=1e4) {
  n <- nrow(X)
  p <- ncol(X)

  mis.ri <- NULL
  k <- 1
  for (i in 1:n) {
    if (any(is.na(X[i,]))) {
      mis.ri[k] <- i
      k <- k + 1
    }
  }

  new.X <- X
  for (j in 1:p) new.X[which(is.na(X[,j])),j] <- mean(X[,j],na.rm=T)

  old.X <- new.X + thresh*2
  j <- 1

  #while (!(all(abs((new.X-old.X)/old.X)<thresh)) & (j<maxit)) {
  while (!(all(abs(new.X-old.X)<thresh)) & (j<maxit)) {
    for (i in mis.ri) {
      if (any(is.na(X[i,]))) {
        mi <- which(is.na(X[i,]))
        mu1 <- apply(as.matrix(new.X[-i, mi]),2,mean)
        mu2 <- apply(as.matrix(new.X[-i,-mi]),2,mean)
        x2  <- new.X[i,-mi]
        S11 <- var(new.X[-i,mi])
        S12 <- var(new.X[-i,mi],new.X[-i,-mi])
        S22 <- var(new.X[-i,-mi])
        B <- S12 %*% solve(S22)
        x1 <- mu1 + B %*% (x2-mu2)
        old.X <- new.X
        new.X[i,mi] <- as.vector(x1)
      }
      #print(Sys.time())
      #print(new.X)
      #Sys.sleep(1)
    } 
    j <- j+1
  } 
  
  print(paste(ifelse(j<maxit,"Converged.","Not Converged."),"Number of Iterations:", j)); cat("\n")
  new.X
}

#X <- EM(dat,thresh=10^(-5),maxit=1000)

em.o2 <- EM(o2,thresh=1e-5,1e4)
em.o4 <- EM(o4,thresh=1e-5,1e4)

mu02 <- matrix(c(1300,120,265,7310,820,45,65,28))
mu04 <- matrix(c(1230,105,275,7360,830,41,75,38))

hot.T2 <- function(Y,mu) {
  mu <- matrix(mu)
  Y <- as.matrix(Y)
  S <- cov(Y)
  n <- nrow(Y)
  p <- ncol(Y)
  y.bar <- apply(Y,2,mean)
  T2 <- n * t(y.bar-mu) %*% solve(S) %*% (y.bar-mu)
  nu <- n-1
  F.stat <- (nu-p+1)/(nu*p) * T2
  p.val <- pf(F.stat,p,nu-p+1,lower=F)
  a <- solve(S) %*% (y.bar-mu)

  info <- matrix(c(T2,n,p,nu,F.stat,p.val),nrow=1)
  colnames(info) <- c("T2","n","p","nu","F.stat","p.val")

  out <- list("info"=info,"discriminant"=a)
  out
}

T2.o2 <- hot.T2(em.o2,mu02)
T2.o4 <- hot.T2(em.o4,mu04)

two.samp.hot.T2 <- function(Y1,Y2) {
  W1 <- (nrow(Y1)-1) * cov(Y1)
  W2 <- (nrow(Y2)-1) * cov(Y2)

  n1 <- nrow(Y1)
  n2 <- nrow(Y2)

  Spl <- (W1+W2)/(n1+n2-2)

  y.bar.1 <- apply(Y1,2,mean)
  y.bar.2 <- apply(Y2,2,mean)
  x <- y.bar.1 - y.bar.2
  T2 <- (n1*n2) / (n1+n2) * t(x) %*% solve(Spl) %*% x

  nu <- n1+n2-2
  p <- ncol(Y1) #REQUIRES: ncol(Y1) = ncol(Y2)
  F.stat <- (nu-p+1)/(nu*p) * T2
  p.val <- pf(F.stat,p,nu-p+1,lower=F)
  a <- solve(Spl) %*% (y.bar.1 - y.bar.2)

  info <- matrix(c(T2,p,nu,F.stat,p.val),nrow=1)
  colnames(info) <- c("T2","p","nu","F.stat","p.val")
  out <- list("info"=info,"discriminant"=a)

  out
}

T2.diff <- two.samp.hot.T2(em.o2,em.o4)

#univariate.tests <- t(apply(matrix(1:8,nrow=1),2,function(x) hot.T2(em.o2[,x],0)))
#colnames(univariate.tests) <- c("T2","n","p","nu","F.stat","p.val")
#univariate.tests

