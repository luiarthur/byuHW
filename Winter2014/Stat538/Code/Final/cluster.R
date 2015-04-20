library(survival)

##################################################
### read in data

dlbcl.data <- as.matrix(read.table("dlbcl_expressiondata.csv",sep=","))
dlbcl.surv <- as.matrix(read.table("dlbcl_survdata.csv",sep=","))
dim(dlbcl.data)
dlbcl.surv
dlbcl.time <- dlbcl.surv[,1]
dlbcl.cens <- dlbcl.surv[,2]

# divide data into training/test sets
set.seed(538)
train.index <- sample(c(1:240),120)
dlbcl.train.data <- dlbcl.data[,train.index]
dlbcl.train.time <- dlbcl.surv[train.index,1]
dlbcl.train.cens <- dlbcl.surv[train.index,2]

dlbcl.test.data <- dlbcl.data[,-train.index]
dlbcl.test.time <- dlbcl.surv[-train.index,1]
dlbcl.test.cens <- dlbcl.surv[-train.index,2]

### hierarchical clustering

dim(t(dlbcl.train.data)) # cluster on rows

# distance matrix
dlbcl.dist <- dist(t(dlbcl.train.data)) # default is euclidean

# linkage
dlbcl.clust <- hclust(dlbcl.dist) # default is complete
plot(dlbcl.clust)  # resultant dendrogram

dlbcl.clust <- hclust(dlbcl.dist,method="ward")
plot(dlbcl.clust)  # resultant dendrogram

#postscript("dlbcl.dendrogram.ps")
#  plot(dlbcl.clust)  
#dev.off()

### HEATMAP: CAN TAKE A LONG TIME #####
#postscript("dlbcl.heatmap.ps",height=10,width=10)
heatmap(t(dlbcl.train.data),Rowv=as.dendrogram(dlbcl.clust),Colv=NA)
#dev.off()
#system("ps2pdf dlbcl.heatmap.ps dlbcl.heatmap.pdf")
#system("rm dlbcl.heatmap.ps")
###########################

cutree(dlbcl.clust,k=2)
dlbcl.2group <- cutree(dlbcl.clust,k=2) - 1
dlbcl.2group
dlbcl.3group <- cutree(dlbcl.clust,k=3)


###### ASSESS PERFORMANCE OF TWO APPROACHES: RISK GROUPS ##########################
# check whether hierarchical clustering can
# successfully distinguish between high-risk / low-risk groups
survdiff(Surv(dlbcl.train.time,dlbcl.train.cens)~c(dlbcl.2group))
coxph(Surv(dlbcl.train.time,dlbcl.train.cens)~as.factor(c(dlbcl.3group)))

### two sample t-test to identify markers whose behavior differs between clusters

group1 <- dlbcl.train.data[,dlbcl.2group==0]
group2 <- dlbcl.train.data[,dlbcl.2group==1]

twosample.ttest <- NA
for(i in 1:nrow(dlbcl.data)) {
  twosample.ttest[i] <- t.test(group1[i,],group2[i,],var.equal=T)$p.value 
}
hist(twosample.ttest)
sum(twosample.ttest < 1e-10)
which(twosample.ttest < 1e-10)
hclust.coef <- which(twosample.ttest < 1e-10)


hclust.res <- coxph(Surv(dlbcl.train.time,dlbcl.train.cens)~t(dlbcl.train.data[hclust.coef,]))
summary(hclust.res)
# need to create a coefficient vector (length=number of markers,
#    all zeros expect for selected markers
hclust.fit.beta <- rep(0,nrow(dlbcl.train.data))
for(i in 1:length(hclust.coef)) {
  hclust.fit.beta[hclust.coef[i]] <- hclust.res$coef[[i]]
}

### check model fit
residuals(hclust.res,type="deviance")
plot(residuals(hclust.res,type="deviance"))



####### assess predictive performance of markers identified from hierarchical clustering
####### on test data set

# risk score for those in training data set
exp(t(hclust.fit.beta)%*%dlbcl.train.data)
# risk score for those in test data set
exp(t(hclust.fit.beta)%*%dlbcl.test.data)
hclust.risk.score.test <- exp(t(hclust.fit.beta)%*%dlbcl.test.data)

# create test set risk groups based on risk score (predicted outcome)
# 0: low-risk group
# 1: high-risk group
risk.group <- ifelse(hclust.risk.score.test < median(hclust.risk.score.test),0,1)

# check whether penalized approach identified markers that
# successfully distinguish between high-risk / low-risk groups
survdiff(Surv(dlbcl.test.time,dlbcl.test.cens)~c(risk.group))

# plot K-M survival functions for high/low risk
g1.time <- dlbcl.test.time[risk.group==0] # low-risk group
g1.cens <- dlbcl.test.cens[risk.group==0]

g2.time <- dlbcl.test.time[risk.group==1] # high-risk group
g2.cens <- dlbcl.test.cens[risk.group==1]

g1.KM <- survfit(Surv(g1.time,g1.cens)~1,type="kaplan-meier")
g2.KM <- survfit(Surv(g2.time,g2.cens)~1,type="kaplan-meier")
summary(g1.KM)
summary(g2.KM)

ymax <- 1
ymin <- 0

g1.KM.time <- c(0,g1.KM$time)
g1.KM.surv <- c(1,g1.KM$surv)
g2.KM.time <- c(0,g2.KM$time)
g2.KM.surv <- c(1,g2.KM$surv)

g1.KM.upper <- c(1,g1.KM$upper)
g1.KM.lower <- c(1,g1.KM$lower)
g2.KM.upper <- c(1,g2.KM$upper)
g2.KM.lower <- c(1,g2.KM$lower)

plot(g1.KM.time,g1.KM.surv,ylim=c(ymin,ymax),type="s",lty=1,
     xlab="Time",
     ylab="Survival",
     main="Unsupervised Hierarchical Clustering: High-Risk vs. Low-Risk",
     cex.main=2,
     cex.axis=1.5,
     cex.lab=1.5,
     xlim=c(0,max(dlbcl.test.time)))
lines(g2.KM.time,g2.KM.surv,lty=1,type="s",col=2)



###### WGCNA clustering approach #####################################################


library(WGCNA)
options(stringsAsFactors = FALSE)
head(dlbcl.train.data)
dlbcl.train.df <- as.data.frame(t(dlbcl.train.data))
sampleTree <- flashClust(dist(dlbcl.train.df), method ="average")
plot(sampleTree)
dlbcl.2group.wgcna = cutree(sampleTree, k=2)

sampleTree <- flashClust(dist(dlbcl.train.df), method ="ward")
plot(sampleTree)
dlbcl.2group.wgcna = cutree(sampleTree, k=2) - 1

survdiff(Surv(dlbcl.train.time,dlbcl.train.cens)~c(dlbcl.2group.wgcna))


