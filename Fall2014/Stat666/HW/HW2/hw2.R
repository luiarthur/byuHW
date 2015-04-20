#$ mail -s "Data" -a "probe.dat" jwardy21@gmail.com

T2.to.F <- function(T2,p,nu) {
  (nu-p+1) / (nu*p) * T2
}

# 12
probe <- read.table("probe.dat",header=T)
n <- nrow(probe)
p <- ncol(probe)

# 12a) Test H0: mu = (30,25,40,25,30)
mu <- matrix(c(30,25,40,25,30),ncol=1)

y.bar <- matrix(apply(probe,2,mean),ncol=1)
y.bar

S <- cov(probe)
S

T2 <- n * t(y.bar-mu) %*% solve(S) %*% (y.bar-mu) 
T2

nu <- n-1
F.stat <- (nu-p+1)/(nu*p) * T2

pf(F.stat,p,nu-p+1,lower=F)

# 12b)

get.t <- function(y,mu) {
  y <- as.vector(y)
  n <- length(y)
  y.bar <- mean(y)
  s <- sd(y)
  sqrt(n) * (y.bar-mu) / s
}

ts <- apply(matrix(1:p,ncol=1),1,function(x) get.t(probe[,x],mu[x,]))
ts <- matrix(ts,nrow=1)
ts

apply(ts,2,function(t) pt(abs(t),n-1,lower=F))

# 16)


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

  out <- matrix(c(T2,n,p,nu,F.stat,p.val),nrow=1)
  colnames(out) <- c("T2","n","p","nu","F.stat","p.val")
  out
}


beet <- as.matrix(read.table("beetles.dat",header=F)[,-1])
ole <- beet[which(beet[,1]==1),][,-1]
car <- beet[which(beet[,1]==2),][,-1]
p <- ncol(ole)

# 16a) Test H0: mu1 = mu2
get.w <- function(Y) {
  n <- nrow(Y)
  (n-1) * cov(Y)
}

W1 <- get.w(ole)
W2 <- get.w(car)

n1 <- nrow(ole)
n2 <- nrow(car)

Spl <- 1/(n1+n2-2) * (W1+W2)

y.bar.1 <- apply(ole,2,mean)
y.bar.2 <- apply(car,2,mean)
x <- y.bar.1 - y.bar.2
T2 <- (n1*n2) / (n1+n2) * t(x) %*% solve(Spl) %*% x

nu <- n1+n2-2
F.stat <- T2.to.F(T2,p,nu)
F.stat

p.val <- pf(F.stat,p,nu-p+1,lower=F)
p.val

# 16b)
a <- solve(Spl) %*% (y.bar.1-y.bar.2)
a

# 20a
goods <- as.matrix(read.table("goods.dat",header=F))
cons <- goods[goods[,1]==1,][,-1]
pros <- goods[goods[,1]==2,][,-1]
p <- ncol(cons)

W1 <- get.w(cons)
W2 <- get.w(pros)

n1 <- nrow(cons)
n2 <- nrow(pros)

Spl <- 1/(n1+n2-2) * (W1+W2)

y.bar.1 <- apply(cons,2,mean)
y.bar.2 <- apply(pros,2,mean)
x <- y.bar.1 - y.bar.2
T2 <- (n1*n2) / (n1+n2) * t(x) %*% solve(Spl) %*% x
T2

nu <- n1+n2-2
F.stat <- T2.to.F(T2,p,nu)
F.stat

p.val <- pf(F.stat,p,nu-p+1,lower=F)
p.val

# 20d
tr <- function(x) sum(diag(x))
nelvan <- function(Y1,Y2) {
  n1 <- nrow(Y1)
  n2 <- nrow(Y2)

  S1 <- cov(Y1)
  S2 <- cov(Y2)
  Se <- S1/n1 + S2/n2

  num <- tr(Se %*% Se) + tr(Se)^2
  den <- 1/(n1-1) * (tr((S1/n1) %*% (S1/n1)) + (tr(S1/n1))^2) +
         1/(n2-1) * (tr((S2/n2) %*% (S2/n2)) + (tr(S2/n2))^2)

  num/den
}

nus <- nelvan(cons,pros)
p <- ncol(cons)

Ts <- function(x1,x2) {
  x1.bar <- matrix(apply(x1,2,mean),ncol=1)
  x2.bar <- matrix(apply(x2,2,mean),ncol=1)

  diff <- x1.bar - x2.bar

  S1 <- cov(x1)
  S2 <- cov(x2)
  n1 <- nrow(x1)
  n2 <- nrow(x2)

  M <- solve(S1/n1 + S2/n2) 
  
  t(diff) %*% M %*% diff
}

Ts <- Ts(cons,pros)

F.stat <- T2.to.F(Ts,p,nus)
F.stat

p.val <- pf(F.stat,p,nu-p+1,lower=F)
p.val

# 21
# 21a
words <- as.matrix(read.table("words.dat",header=T))
D <- words[,5:6]
n <- nrow(D)
p <- ncol(D)
mu <- matrix(0,2,1)

y.bar <- matrix(apply(D,2,mean),ncol=1)
y.bar

S <- cov(D)
S

T2 <- n * t(y.bar-mu) %*% solve(S) %*% (y.bar-mu) 
T2

nu <- n-1
F.stat <- T2.to.F(T2,p,nu)
F.stat

pf(F.stat,p,nu-p+1,lower=F)

# 21b
Spl <- cov(D)
a <- solve(Spl) %*% matrix(apply(D,2,mean))



# 21b
ts <- apply(matrix(1:p,ncol=1),1,function(x) get.t(D[,x],0))
ts <- matrix(ts,nrow=1)
ts

apply(ts,2,function(t) pt(abs(t),n-1,lower=F))

# 6b
a <- matrix(c(3,-2,0,0))
mu <- matrix(c(3,11,-4,1))
Sig <- matrix(c(11,3,3,1,
              3,10,2,-5,
              3,2,10,-6,
              1,-5,-6,15),4,4,byrow=T)

t(a) %*% mu
t(a) %*% Sig %*% a

mu.y <- -13
mu.x <- matrix(c(-4,1))
Sxx <- Sig[3:4,3:4]
Syx <- matrix(c(5,13),nrow=1)
x <- matrix(c(1,-2))
EyGz <- mu.y + Syx %*% solve(Sxx) %*% (x-mu.x)
EyGz

Syy <- 103
CyGz <- Syy - Syx %*% solve(Sxx) %*% t(Syx)
CyGz


# 6c
det(S)

det(diag(diag(S)))

