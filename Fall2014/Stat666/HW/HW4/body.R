# I have questions at #5,#6,#7
#options("width"=150)

lam.to.F <- function(lam,p,vh,ve) { # Approximation
  t <- sqrt( (p^2*vh^2-4) / (p^2+vh^2-5) )
  L <- lam^(1/t)
  w <- ve + vh - (p+vh+1)/2

  df <- c(p*vh, w*t-(p*vh-2)/2)

  F.stat <- (1-L) / L * df[2]/df[1]
  p <- pf(F.stat,df[1],df[2],lower.tail=F)
  out <- c(F.stat,df,p)
  names(out) <- c("F.stat","df1","df2","p.val")

  out
}

V.2.F <- function(V,s,N,m) { # Approximation
  F.stat <- (2*N + s + 1) / (2*m + s + 1) * V/(s-V)
  df <- s * (c(m,N) + s + 1)
  p <- pf(F.stat,df[1],df[2],lower.tail=F)

  out <- c(F.stat,df,p)
  names(out) <- c("F.stat","df1","df2","p")
  
  out
}

#5:
dat <- as.matrix(read.table("bodyfat.txt",header=F))
colnames(dat) <- c(paste0("y",1:2), paste0("x",1:13))

#y1  = Density determined from underwater weighing
#y2  = Percent body fat
#x1  = Age (years)
#x2  = Weight (lbs)
#x3  = Height (inches)
#x4  = Neck circumference (cm)
#x5  = Chest circumference (cm)
#x6  = Abdomen 2 circumference (cm)
#x7  = Hip circumference (cm)
#x8  = Thigh circumference (cm)
#x9  = Knee circumference (cm)
#x10 = Ankle circumference (cm)
#x11 = Biceps (extended) circumference (cm)
#x12 = Forearm circumference (cm)
#x13 = Wrist circumference (cm)

Y <- dat[,1:2]
X <- dat[,-(1:2)]

n <- nrow(Y)
p <- ncol(Y)
q <- ncol(X)
r <- q+1

X <- cbind(1,X)
B <- solve(t(X)%*%X) %*% t(X)%*%Y

E <- t(Y-X%*%B) %*% (Y-X%*%B)
S <- E / (n-r)
y.bar <- apply(Y,2,mean)
H <- t(Y) %*% Y - n*y.bar%*%t(y.bar)
lam <- det(E) / det(E+H)


s <- min(p,q)
l <- eigen(solve(E,H))$values[1:s]
#lam <- prod(1/(1+l))
V <- sum(1/(1+l))

F.L <- lam.to.F(lam,p=p,vh=q-1,ve=n-r)
F.V <- V.2.F(V,s=min(p,q),N=(n-q-p-2)/2,m=(abs(p-1)-1)/2)


Syy <- var(Y)
Sxx <- var(X[,-1])
Syx <- cov(Y,X[,-1])
Sxy <- t(Syx)
A <- solve(Syy) %*% Syx %*% solve(Sxx) %*% Sxy
r2 <- eigen(A)$values[1:s] # Squared Canonical Correlation in SAS
                           # the sqrt of r2 is the canonical corr in SAS
R2 <- prod(r2) # Should be the same as det(A)
#R2 <- det(A) # Should be the same as product of the eigen values of A

L.m <- apply(matrix(1:s),1,function(m) prod(1-r2[m:length(r2)]))
F.m <- t(apply(matrix(1:s),1,function(m) lam.to.F(L.m[m],p=q-m+1,vh=p-m+1,ve=n-m-p)))

r2.l <- r2/(1-r2) # This is the "eigenvalues" of E^(-1)H in SAS = l-1. Why?
prop.r2.l <- r2.l / sum(r2.l) # 98% => essentially one dimension 

#6:
back.sel <- function(Y,X,a=.05) {
  n <- nrow(Y)
  p <- ncol(Y)
  X <- cbind(1,X)
  YTY <- t(Y) %*% Y

  get.L <- function(ve) {# note: p=p, vh=1, a=a=.05
    F.stat <- qf(1-a,p,ve-p+1)
    L <- 1 / (F.stat*p/(ve-p+1)+1)
    L
  }

  get.L.xi.rest <- function(i,x,E.full) {
    x <- x[,-i]
    b <- solve(t(x)%*%x, t(x)%*%Y)
    E.red <- YTY - t(x%*%b) %*% Y
    L <- det(E.full)/det(E.red)
    L
  }

  L.end <- get.L(n-r)
  L <- L.end + 1

  while (L > L.end) {
    r <- ncol(X)
    B <- solve(t(X)%*%X, t(X)%*%Y)
    E.full <- YTY - t(X%*%B) %*% Y
    #L.s <- apply(matrix(2:r),1,function(i) get.L.xi.rest(i,X,E.full))
    L.s <- apply(matrix(1:r),1,function(i) get.L.xi.rest(i,X,E.full))
    i <- which.max(L.s)
    L <- L.s[i]
    L.end <- get.L(n-r)
    if (L > L.end) {
      #X <- X[,-(i+1)]
      X <- X[,-i]
    } 
  }

  X
}

new.X <- back.sel(Y,X[,-1])
#head(new.X)

# Question: Do I need to consider removing the intercept, or do I always
#           keep the intercept?

#7:
Sxx <- var(new.X[,-1])
Syx <- cov(Y,new.X[,-1])
Sxy <- t(Syx)
A <- solve(Syy) %*% Syx %*% solve(Sxx) %*% Sxy
r2 <- eigen(A)$values[1:min(ncol(new.X),ncol(Y))] # Squared Canonical Correlation in SAS
                           # the sqrt of r2 is the canonical corr in SAS
R2 <- prod(r2) # Should be the same as det(A)
#R2 <- det(A) # Should be the same as product of the eigen values of A

L.m <- apply(matrix(1:s),1,function(m) prod(1-r2[m:length(r2)]))
F.m <- t(apply(matrix(1:s),1,function(m) lam.to.F(L.m[m],p=q-m+1,vh=p-m+1,ve=n-m-p)))

r2.l <- r2 / (1-r2)
r2.l / sum(r2.l) # 99% in eig[1]

B.new <- solve(t(new.X) %*% new.X, t(new.X)%*%Y) 



Syy <- var(Y)
Syx <- var(Y,new.X[,-1])
Sxy <- var(new.X[,-1],Y)
Sxx <- var(new.X[,-1])

A <- solve(Syy,Syx) %*% solve(Sxx,Sxy)
B <- solve(Sxx,Sxy) %*% solve(Syy,Syx)

a <- eigen(A)$vector
b <- eigen(B)$vector

(c <- diag(sqrt(diag(Syy))) %*% a) 
(d <- diag(sqrt(diag(Sxx))) %*% b) 

rownames(c) <- c("y1","y2")
rownames(d) <- c("x2","x6","x11","x13")
colnames(c) <- c("c1","c2")
colnames(d) <- c("d1","d2","d3","d4")

c <- c[,1:min(ncol(c),ncol(d))]
d <- d[,1:min(ncol(c),ncol(d))]
# Why is my answer different???
