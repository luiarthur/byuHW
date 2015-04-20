library(survival)
source("../../Midterm/PCNSL/R/km.perc.R") #se.pcnt
##################################################
# Clean Data:
  # Expression Data:
  exprss <- read.csv("../Data/BladderCancer_expression.csv")
  # Make the first column the column name
  rownames(exprss) <- exprss[,1]     # These are the gene names
  exprss <- exprss[-nrow(exprss),-1] # Remove the gene names and column of NA's

  # Clinic Data:
  clinic <- read.csv("../Data/BladderCancer_clinical.csv")
  clinic <- clinic[,c("overall.survival","survivalMonth")] # remove patient ID's
  colnames(clinic) <- c("cens","time") 
  clinic$cens <- clinic$cens - 1 # 1 = censored, 2 = death
                                 # Want: 0 = cens, 1 = death

  # No Covariates:
    bigD <- cbind(clinic,t(exprss))
    N <- nrow(bigD)
    P <- ncol(bigD) - 2
    censor.rate <- mean(bigD$cens == 0)
    #pdf("latex/raw/cens.pdf")
    #  plot(bigD[,1:2],col='red')
    #  legend("center",legend=c("Censored = 0","Death = 1") )
    #dev.off()  

    time.cens <- bigD[,2][which(bigD[,1]==0)]
    time.died <- bigD[,2][which(bigD[,1]==1)]
    time.summary <- rbind( summary(time.cens),
                           summary(time.died), 
                           summary(bigD[,2]) )
    rownames(time.summary) <- c("Censored","Died","Overall")
    #library(xtable)
    #sink("latex/raw/timeSummary.tex")
    #  xtable(time.summary)
    #sink()

    #pdf("latex/raw/times.pdf")
    #  par(mfrow=c(3,1),mar=rep(6,4))
    #  hist(bigD[,2],main="Histogram of Time",xlab="Time (Days)",col='yellow',br=50)
    #  hist(time.cens,main="Histogram of Censored Time",
    #       xlab="Censored Time (Days)",col='yellow',br=50)
    #  hist(time.died,main="Histogram of Death Time",
    #       xlab="Death Time (Days)",col='yellow',br=50)
    #  par(mfrow=c(1,1))
    #dev.off()  

# divide data into training/test sets
  set.seed(1)
  #train.ind <- sample(1:N,round(N/2))
  train.ind <- sample(1:N,130)
  train.set <- bigD[ train.ind,]
  test.set  <- bigD[-train.ind,]
  train.exprss <- train.set[,-c(1:2)]
  test.exprss  <-  test.set[,-c(1:2)]
  
### hierarchical clustering

dim(train.exprss) # cluster on rows

# distance matrix
dist <- dist(train.exprss,method="maximum") # default is euclidean # Devyn did maximum

# linkage
clust <- hclust(dist,method="ward") # Devyn used complete
#pdf("latex/raw/ward.pdf")
#  plot(clust)  # resultant dendrogram
#dev.off()

### HEATMAP: CAN TAKE A LONG TIME #####
#pdf("latex/raw/heatmap.pdf")
#heatmap(as.matrix(train.exprss),Rowv=as.dendrogram(clust),Colv=NA)
#dev.off()
# No idea what we're looking for in this plot...

###########################
#####???
cutree(clust,k=2)
clust.2group <- cutree(clust,k=2) - 1
clust.3group <- cutree(clust,k=3) - 1

###### ASSESS PERFORMANCE OF TWO APPROACHES: RISK GROUPS ##########################
# check whether hierarchical clustering can
# successfully distinguish between high-risk / low-risk groups
x <- survdiff(Surv(train.set$time,train.set$cens)~c(clust.2group))
#sink("latex/raw/clust2.pdf")
#  x
#sink()

x <- coxph(Surv(train.set$time,train.set$cens)~as.factor(c(clust.3group)))
#sink("latex/raw/clust3.pdf")
#  xtable(x)
#  cat("\\tiny Likelihood ratio test=4.05  on 2 df, p=0.132  n= 130, number of events= 49","\n")
#sink()  

### two sample t-test to identify markers whose behavior differs between clusters

group1 <- train.exprss[clust.2group==0,]
group2 <- train.exprss[clust.2group==1,]

####???express or bigD
twosample.ttest <- NA

for(i in 1:N) {
  twosample.ttest[i] <- t.test(group1[,i],group2[,i],var.equal=T)$p.value 
}
#pdf("latex/raw/histClust.pdf")
#  hist(twosample.ttest)
#dev.off()  
sum(twosample.ttest < 1e-5)
which(twosample.ttest < 1e-5)
hclust.coef <- which(twosample.ttest < 1e-5)

genes <- as.matrix(train.exprss[,hclust.coef])
colnames(genes) 
hclust.res <- coxph(Surv(train.set$time,train.set$cens) ~ genes)
summary(hclust.res)
#sink("latex/raw/clustCoef.tex")
#  xtable(hclust.res)
#sink()

# need to create a coefficient vector (length=number of markers,
#    all zeros expect for selected markers
hclust.fit.beta <- rep(0,ncol(train.exprss))
for(i in 1:length(hclust.coef)) {
  hclust.fit.beta[hclust.coef[i]] <- hclust.res$coef[[i]]
}

### check model fit
residuals(hclust.res,type="deviance")
#pdf("latex/raw/devResTree.pdf")
#  plot(residuals(hclust.res,type="deviance"))
#dev.off()

# assess predictive performance of markers identified from hierarchical clustering
# on test data set

# risk score for those in training data set
exp(as.matrix(train.exprss)%*%hclust.fit.beta)
# risk score for those in test data set
hclust.risk.score.test <- exp(as.matrix(test.exprss)%*%hclust.fit.beta)

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
risk.group <- ifelse(hclust.risk.score.test < median(hclust.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
library(rms)

risk <- ifelse(risk.group==1,"High","Low")
#pdf("latex/raw/clustKM.pdf")
#  survplot(km <- survfit(Surv(test.set$time,test.set$cens)~risk))
#  title(main="Unsupervised Hierarchical Clustering: High-Risk vs. Low-Risk")
#dev.off()

lo <- which(risk=="Low")
hi <- which(risk=="High")
km.hi <- survfit(Surv(test.set$time[hi],test.set$cens[hi])~risk[hi])
km.lo <- survfit(Surv(test.set$time[lo],test.set$cens[lo])~risk[lo]) 
se.pcnt(.5,km.hi)
se.pcnt(.5,km.lo)
#km

##############################################

library(survivalROC)

time.seq <- seq(0,140,by=.5)
clust.AUC.res <- NA
for(i in 1:length(time.seq)) {
clust.AUC.res[i] <- survivalROC(Stime=test.set$time,status=test.set$cens,marker=hclust.risk.score.test,predict.time=time.seq[i],method="KM")$AUC
}
#write.table(clust.AUC.res,"aucData/clust.txt",quote=F,col.names=F,row.names=F)

pdf("latex/raw/clustAUC.pdf")
plot(time.seq,clust.AUC.res,type="l",ylim=c(0,1),main="Time-Dependent ROC",
     ylab="AUC",xlab="Time",col="pink",lwd=3)
dev.off()


# Just Run This Chunk of Code: ##############################################
#time.seq <- seq(0,140,by=.5)
#clust.auc <- as.matrix(read.table("aucData/clust.txt"))
#margi.auc <- as.matrix(read.table("aucData/marginal.txt"))
#pen.auc   <- as.matrix(read.table("aucData/pen.txt"))
#rf.auc    <- as.matrix(read.table("aucData/rf.txt"))
#pca.auc   <- as.matrix(read.table("aucData/pca.txt"))
#
#final.plot <- function(cx = 1) {
#  par(mfrow=c(2,1))
#  plot(time.seq,clust.auc,type="l",ylim=c(0,1),main="Time-Dependent ROC",
#       ylab="AUC",xlab="Time",col="pink",lwd=3,xlim=c(0,140))
#  lines(time.seq,margi.auc,type="l",col="blue",lwd=3)
#  lines(time.seq,pen.auc,type="l",col="red",lwd=3)
#  lines(time.seq,rf.auc,type="l",col="green",lwd=3)
#  lines(time.seq,pca.auc,type="l",col="orange",lwd=3)
#  legend("topright",legend=c("Lasso","FDR","H-Cluster","Random Forest","PCA"),
#                    col=c("red","blue","pink","green","orange"),
#                    title="Method Used",lwd=3,cex=cx)
#  hist(bigD$time,breaks=50,xlab="Time",main="Histogram of Time",col="yellow")
#  par(mfrow=c(1,1))
#}
#pdf("latex/raw/finalPlot.pdf"); final.plot(.6); dev.off()

