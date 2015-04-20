# 3.11
Y <- read.table("data",header=T)[,-1]
Y <- as.matrix(Y)

# 3.11a
S <- cov(Y)
dS <- det(S)
dS

#3.11b
tr <- function(M) sum(diag(M)) # where M is a square matrix
trS <- tr(S)
trS

#3.15
w <- matrix(c(-2,3,1),ncol=1)
z <- matrix(c(3,-1,2),ncol=1)

z.vec <- matrix(c(apply(Y,1,function(x) x %*% z)))
colnames(z.vec) <- "z"
z.vec

w.vec <- matrix(c(apply(Y,1,function(x) x %*% w)))
colnames(w.vec) <- "w"
w.vec

# a)
s.xy <- function(x,y) {
  m.x <- mean(x)
  m.y <- mean(y)
  n <- length(x)

  sum((x-m.x)*(y-mean(y))) / (n-1)
}

s <- function(x) {
  m.x <- mean(x)
  n <- length(x)

  sqrt(sum((x-m.x)^2/(n-1)))
}

r <- function(x,y) {
  s.xy(x,y) / (s(x) * s(y))
}

zw <- cbind(z.vec,w.vec)
r(zw[,1],zw[,2])

# b)
a <- z 
b <- w

t(a) %*% S %*% b / sqrt((t(a) %*% S %*% a) %*% (t(b) %*% S %*% b))


# 3.17)
# a)
a <- matrix(c(1,1,1),ncol=1)
b <- matrix(c(2,-3,2),ncol=1)
c <- matrix(c(-1,-2,-3),ncol=1)

A <- rbind(t(a),t(b),t(c))
z.bar <- A %*% apply(Y,2,mean)
z.bar

S.z <- A %*% S %*% t(A)
S.z

# b)
D.s <- diag(sqrt(diag(S.z)))
R.z <- solve(D.s) %*% S.z %*% solve(D.s)
R.z


# 3.22)
glu <- read.table("glucose.dat",header=T)

Y <- glu[,1:3]
X <- glu[,4:6]

y.bar <- apply(Y,2,mean)
x.bar <- apply(X,2,mean)

mean.vec <- matrix(c(y.bar,x.bar),ncol=1)
mean.vec

S <- cov(glu)
S


# 4.1
S1 <- matrix(c(14,8,3,8,5,2,3,2,1),nrow=3)
S2 <- matrix(c(6,6,1,6,8,2,1,2,1),nrow=3)

det(S1)
det(S2)

D1 <- diag(sqrt(diag(S1)))
D2 <- diag(sqrt(diag(S2)))

R1 <- solve(D1) %*% S1 %*% solve(D1)
R2 <- solve(D2) %*% S2 %*% solve(D2)

D1
D2

R1
R2

# 4.11
# a)
S <- matrix(c(6,1,-2,1,13,4,-2,4,4),nrow=3)
M <- chol(S)
M

solve(t(M))

# b)
D.sqrt <- diag(sqrt(eigen(S)$values))
C <- eigen(S)$vectors
S.sqrt <- C %*% D.sqrt %*% t(C)
solve(S.sqrt)


# 4.17)
mu <- matrix(c(-3,2,4,-3,5),ncol=1)
S.yx <- matrix(c(15,0,3,8,6,-2),nrow=2,byrow=T)
S.xx <- matrix(c(50,8,5,8,4,0,5,0,1),ncol=3)
S.yy <- matrix(c(14,-8,-8,18),ncol=2)

#a
S.yx %*% solve(S.xx) 

#b
S.yy - S.yx %*% solve(S.xx) %*% t(S.yx)


# 11
Y <- read.table("glucose.dat",header=T)[1:10,1:5]
y <- read.table("glucose.dat",header=T)[11,1:5]
Y.new <- rbind(Y,y)

An.Inv <- solve(cov(Y))

A.np1.1 <- solve(cov(Y.new))

update.An.Inv <- function(An.Inv,Y,y,n) {
  Y <- as.matrix(Y)

  sherman <- function(Binv,c) {
    # (B + cc')^-1:
    Binv - 
    Binv %*% c %*% t(c) %*% Binv / 
    as.numeric(1 + t(c) %*% Binv  %*% c)
  }
  
  y.bar <- matrix(apply(as.matrix(Y),2,mean),ncol=1)
  y <- Y[nrow(Y),]
  p <- nrow(An.Inv)
  Binv <- n/(n-1) * An.Inv 
  c <- sqrt(1+n) / n * 
  (y.bar - y)

  sherman(Binv,c)
}

A.np1.2 <- update.An.Inv(An.Inv,Y.new,n=10)

A.np1.1
A.np1.2
