# SectVC.R
source("multT2.R")
lrbind <- function(L) { #rbinds a list of matrices
  k <- length(L)
  p <- ncol(L[[1]])
  out <- matrix(0,0,p)
  for (i in 1:k) {
    out <- rbind(out,L[[i]])
  }
  out
}

X <- read.table("collins.txt",header=T)
write.table(X,quote=F,col.names=F,row.names=F,file="cleanData.txt")

rownames(X) <- X[,1]
Y <- X[,-1]

Z <- scale(Y[,1:18])
Z <- data.matrix(Z)
D <- dist(Z)
wardlink <- hclust(D,method="ward.D")
plot(wardlink,labels=rownames(Y),cex=.0001)

grps <- as.list(3:7)
grps <- lapply(grps,function(x) cutree(wardlink,k=x))

centers <- as.list(3:7)
z <- as.list(3:7)
for (i in 1:5) {#3:7
  z[[i]] <- as.list(1:(i+2))
  centers[[i]] <- matrix(0,i+2,ncol(Z))
  for (j in 1:(i+2)){
    # e.g. z[[1]][[2]] <- for k=3, get the 2nd cluster
    z[[i]][[j]] <- Z[which(grps[[i]]==j),]
    centers[[i]][j,] <- apply(z[[i]][[j]],2,mean)
  }
}

km <- as.list(3:7)
for (i in 3:7) {
  km[[i-2]] <- kmeans(Z,centers[[i-2]])
}
km.grps <- lapply(km,function(x) x$cluster) #k=3:7

z.new <- as.list(1:5)
for (i in 3:7) {
  z.new[[i-2]] <- as.list(1:i)
  for (j in 1:i){
    z.new[[i-2]][[j]] <- Z[which(km.grps[[i-2]]==j),] #data for k=3, 1st cluster
  }
}

#ybars <- zbar[[3]]
#ys <- z[[3]]
my.manova <- function(ys) {# a list of matrices
  k <- length(ys)
  n <- as.numeric(lapply(ys,nrow))
  p <- ncol(ys[[1]])
  N <- sum(n)
  ybar. <- lapply(ys,function(y) as.matrix(apply(y,2,mean)))
  y <- lrbind(ys)
  ybar.. <- apply(y,2,mean)

  H <- matrix(0,p,p)
  for (i in 1:k) {
    A <- ybar.[[i]]-ybar..
    H <- n[i] * A%*%t(A) + H
  }
  
  E <- matrix(0,p,p)
  for (i in 1:k){
    for (j in 1:(n[i])){
      A <- ys[[i]][j,] - ybar.[[i]]
      E <- A%*%t(A) + E
    }
  }
  
  lam <- det(E) / det(E+H)
  vH <- k-1
  vE <- N-k
  
  t <- sqrt(((p*vH)^2-4)/(p^2+vH^2-5))
  w <- vE+vH-.5*(p+vH+1)
  df1 <- p*vH
  df2 <- w*t-.5*(p*vH-2)
  lam.th.rt <- lam^(1/t)
  F.stat <- (1-lam.th.rt)*df2/(lam.th.rt*df1)
  #see p.185 of 666 (373)
  out <- c(F.stat,df1,df2,pf(F.stat,df1,df2,lower.tail=F))
  names(out) <- c("F.stat","df1","df2","p.val")
  out
}

manova.result <- matrix(0,5,4)
for (i in 3:7) {
  manova.result[i-2,] <- my.manova(z.new[[i-2]])
}

colnames(manova.result) <- c("F.stat","df1","df2","p.val")
rownames(manova.result) <- paste0("k=",3:7)
manova.result


# Test predict.kmean
library(clue)
cl_predict(km[[1]],Z[999:1000,])

