# Code:

source("multT2.R")
source("a.manova.R")
source("lrbind.R")
source("rapply.R")
source("countdown.R")
source("a.image.R")

X <- read.table("collins.txt",header=T)
write.table(X,quote=F,col.names=F,row.names=F,file="cleanData.txt")

rownames(X) <- X[,1]
Y <- X[,-1]
Z <- scale(Y[,1:18])

get.clusters <- function(Z,k=3,meth="ward.D") { #z=standardized data, k=num of clusters
  D <- dist(Z)
  link <- hclust(D,method=meth)
  grp <- cutree(link,k=k)

  z <- as.list(1:k)
  center <- matrix(0,k,ncol(Z))
  for (j in 1:k) {
    z[[j]] <- Z[which(grp==j),]
    center[j,] <- apply(z[[j]],2,mean)
  }
  
  km <- kmeans(Z,center)
  km.grp <- km$cluster

  z.new <- as.list(1:k)
  for (j in 1:k) {
    z.new[[j]] <- Z[which(km.grp==j),]
  }

  out <- list("km"=km,"z"=z.new) # km = kmeans object, z = new data list
  out
}

ks <- 3:7 # number of clusters
clusters <- lapply(as.list(ks),function(x) get.clusters(Z,x,"ward.D2"))
#clusters <- lapply(as.list(ks),function(x) get.clusters(Z,x,"complete"))
clus.z <- lapply(clusters,function(x) x$z)
clus.km <- lapply(clusters,function(x) x$km)

manova.result <- rapply(clus.z,a.manova)

# Changes:
  temp <- a.manova(clus.z[[1]])
  E.in.H <- solve(temp$E,temp$H)
  eig.vec <- Re(eigen(E.in.H)$vec[,1:2])

  x <- Z %*% eig.vec
  library(xtable)
  rownames(eig.vec) <- colnames(Z)
  sink("EH.tex")
    xtable(eig.vec)
  sink()


manova.stat <- rapply(manova.result,function(x) x$stats)
colnames(manova.stat) <- c("F.stat","df1","df2","p.val")
rownames(manova.stat) <- paste0("k=",3:7)
manova.stat


#V2:
#Changes:
  par(mfrow=c(1,1))
pdf("disc.pdf")
  plot(x, type="n", xlab="Discriminant Func. Score 1", 
                    ylab="Discriminant Func. Score 2",main="K-means")
  text(x, labels=clus.km[[k]]$cluster, col=clus.km[[k]]$cluster,cex=.7,lwd=9)
  #text(x, labels=labels(Z)[[1]],col=clus.km[[k]]$cluster,cex=.6)
  #points(cen[,1:2], pch=3, cex=3,col="blue",lwd=3)
dev.off()


pca <- princomp(Z)
px <- predict(pca)
k <- 1
cen <- predict(pca,clus.km[[k]]$centers)
pdf("pca.pdf")
  plot(px[,1:2], type="n", xlab="PC1", ylab="PC2",main="K-means")
  text(px[,1:2], labels=clus.km[[k]]$cluster, col=clus.km[[k]]$cluster,cex=.7,lwd=9)
  #text(px[,1:2], labels=labels(Z)[[1]],col=clus.km[[k]]$cluster,cex=.6)
  points(cen[,1:2], pch=3, cex=3,col="blue",lwd=3)
dev.off()


library(clue) # for cl_predict
clus.cv <- function(Z,k=3) { # Z is your standardized data, k is # of clusters
  all.km <- get.clusters(Z,k)$km
  n <- nrow(Z)
  out <- NULL

  one.it <- function(i) {
    Z.i <- Z[-i,]
    one.out.km <- get.clusters(Z.i,k)$km
    cl_predict(one.out.km,matrix(Z[i,],1)) == all.km$cluster[i]
  }

  #V1:
  #apply(matrix(1:n),1,one.it)

  #V2:
  out <- NULL
  for (i in 1:n) {
    old.time <- Sys.time()
      out[i] <- one.it(i)
    count.down(old.time,i,n)
  }
  out
}

library(doMC)
registerDoMC(strtoi(system("nproc",intern=T))/2)
error.rate <- foreach(i=1:length(ks),.combine=cbind) %dopar% (1-mean(clus.cv(Z,ks[i])) )
colnames(error.rate) <- paste0("k=",ks)
rownames(error.rate) <- "Error.Rate"

#  error.rate
#  k=3   k=4   k=5   k=6   k=7 
#  0.008 0.234 0.245 0.307 0.330 

#1: Press
#2: Non-press
#3: Biography
#4: Scholarship
#5: Fiction
gen <- Y$Genre
supGen <- ifelse(gen %in% 1:3,1, ifelse(gen %in% 4:6,2, ifelse(gen == 7,3, ifelse(gen %in% 8:9,4,5))))

sup.Z <- lapply(as.list(unique(supGen)),function(x) Z[which(supGen==x),])
sup.manova <- a.manova(sup.Z)
sup.manova.stat <- sup.manova$stats
comp.sup.k5 <- rbind(manova.stat[3,],sup.manova.stat)
rownames(comp.sup.k5) <- c("Natural, k=5:","Super Genres:")
sink("compSupK5.tex")
  xtable(comp.sup.k5)
sink()  


supGen.k3.cluster <- table(supGen,clus.km[[1]]$cluster)
supGen.prop <- t(apply(supGen.k3.cluster,1,function(x) x/sum(x)))
rownames(supGen.prop) <- c("Press","Non-press","Biography","Scholarship","Fiction")
colnames(supGen.prop) <- paste0("Cluster",1:3)
supGen.prop
a.image(supGen.prop,cex.axis=.8,lasx=0,col=rev(heat.colors(20)))

gen.k3.cluster <- table(gen,clus.km[[1]]$cluster)
gen.prop <- t(apply(gen.k3.cluster,1,function(x) x/sum(x)))
colnames(gen.prop) <- colnames(supGen.prop)
rownames(gen.prop) <- c("Press:\nReporting","Press:\nEditorial","Press:\nReviews",
                        "Religion","Skills\n& Hobbies","Popular\nLore","Biography",
                        "Official\nComm.","Learned","General\nFiction",
                        "Mystery","Science\nFiction","Adventure","Romance","Humor")
gen.prop
a.image(gen.prop,cex.axis=.8,col=rev(heat.colors(20)),lasx=0)


supGen.k5.cluster <- table(supGen,clus.km[[3]]$cluster)
supGen.k5.prop <- t(apply(supGen.k5.cluster,1,function(x) x/sum(x)))
rownames(supGen.k5.prop) <- c("Press","Non-press","Biography","Scholarship","Fiction")
colnames(supGen.k5.prop) <- paste0("Cluster",1:5)
supGen.k5.prop
pdf("supK5prop.pdf")
a.image(supGen.k5.prop,cex.axis=.8,lasx=0,col=rev(heat.colors(20)))
dev.off()

sink("supGenk5.tex")
  xtable(supGen.k5.prop)
sink()

gen.k5.cluster <- table(gen,clus.km[[3]]$cluster)
gen.k5.prop <- t(apply(gen.k5.cluster,1,function(x) x/sum(x)))
colnames(gen.k5.prop) <- paste0("Cluster",1:5)
rownames(gen.k5.prop) <- c("Press:\nReporting","Press:\nEditorial","Press:\nReviews",
                        "Religion","Skills\n& Hobbies","Popular\nLore","Biography",
                        "Official\nComm.","Learned","General\nFiction",
                        "Mystery","Science\nFiction","Adventure","Romance","Humor")
gen.k5.prop
pdf("../latex/join.pdf",width=14)
par(mfrow=c(1,2))
  a.image(supGen.k5.prop,cex.axis=.8,lasx=0,col=rev(heat.colors(20)))
  a.image(gen.k5.prop,cex.axis=.8,col=rev(heat.colors(20)),lasx=0)
par(mfrow=c(1,1))
dev.off()


cv.supGen <- function(Z) {
  euclid.dist <- function(x,y) sum((x-y)^2)

  one.out.centroid <- function(i) {
    z <- as.list(unique(supGen))
    z <- lapply(z,function(x) Z[which(supGen[-i]==x),])
    centroid <- rapply(z,function(x) apply(x,2,mean))
    centroid
  }

  n <- nrow(Z) 
  clust <- NULL
  for (i in 1:n) {
    centroid <- one.out.centroid(i)
    clust[i] <- which.min(apply(centroid,1,function(x) euclid.dist(x,Z[i,])))
  }
  1-mean(clust==supGen)
}

error.rate.supGen <- cv.supGen(Z)
err.natural.v.supGen <- c(error.rate.supGen,error.rate[3])
err.natural.v.supGen <- matrix(err.natural.v.supGen,1)
colnames(err.natural.v.supGen) <- c("Super.Genres","Natural,k=5")
rownames(err.natural.v.supGen) <- c("Error Rate")
err.natural.v.supGen
#  err.natural.v.supGen
#  Super.Genres  Natural,k=5 
#         0.363        0.245

sink("errRate.tex")
  xtable(error.rate)
sink()  

sink("errVS.tex")
  xtable(err.natural.v.supGen)
sink()


