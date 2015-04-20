
# Just Run This Chunk of Code: ##############################################
time.seq <- seq(0,140,by=.5)
clust.auc <- as.matrix(read.table("aucData/clust.txt"))
margi.auc <- as.matrix(read.table("aucData/marginal.txt"))
pen.auc   <- as.matrix(read.table("aucData/pen.txt"))
rf.auc    <- as.matrix(read.table("aucData/rf.txt"))
pca.auc   <- as.matrix(read.table("aucData/pca.txt"))

final.plot <- function(cx = 1) {
  par(mfrow=c(2,1))
  plot(time.seq,clust.auc,type="l",ylim=c(0,1),main="Time-Dependent ROC",
       ylab="AUC",xlab="Time",col="pink",lwd=3,xlim=c(0,140))
  lines(time.seq,margi.auc,type="l",col="blue",lwd=3)
  lines(time.seq,pen.auc,type="l",col="red",lwd=3)
  lines(time.seq,rf.auc,type="l",col="green",lwd=3)
  lines(time.seq,pca.auc,type="l",col="orange",lwd=3)
  legend("topright",legend=c("Lasso","FDR","H-Cluster","Random Forest","PCA"),
                    col=c("red","blue","pink","green","orange"),
                    title="Method Used",lwd=3,cex=cx)
  hist(bigD$time,breaks=50,xlab="Time",main="Histogram of Time",col="yellow")
  par(mfrow=c(1,1))
}
pdf("latex/raw/finalPlot.pdf"); final.plot(.5); dev.off()


