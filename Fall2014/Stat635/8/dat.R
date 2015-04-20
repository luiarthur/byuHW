dat <- as.matrix(read.csv("primate.dat",header=T))
n <- nrow(dat)
new <- matrix(0,n*5,4)

colnames(new) <- c("Monkey","Treatment","Week","Score")
new[,1] <- rep(dat[,1],each=5)
new[,2] <- rep(dat[,2],each=5)
new[,3] <- rep(c(2,4,8,12,16),n)
new[,4] <- c(t(dat[,3:7]))

write.table(new,file="prim2.txt",quote=F,col.names=T,row.names=F)
