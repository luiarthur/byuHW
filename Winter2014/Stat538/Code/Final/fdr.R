rm(list=ls())
library(survival)

# generate data from two distinct populations
reject.res <- p.value <- NA
for (i in 1:10000) {
  x <- rnorm(100,4,3)
  y <- rnorm(100,5,3)
  p.value[i] <- t.test(x,y,"two.sided")$p.value
  reject.res[i] <- ifelse(p.value[i] < .05,1,0)
}
mean(reject.res)
hist(p.value,breaks=seq(0,1,by=.05))


# generate data from two identical populations
reject.res <- p.value <- NA
for (i in 1:10000) {
  x <- rnorm(1000,5,3)
  y <- rnorm(1000,5,3)
  p.value[i] <- t.test(x,y,"two.sided")$p.value
  reject.res[i] <- ifelse(p.value[i] < .05,1,0)
}
mean(reject.res) # Probability of Type I Error
hist(p.value,breaks=seq(0,1,by=.05))

# in case where populations are the same, will have lots of false positives
# in case where populations are not, will still have some false positives


##################################################
### read in data

dlbcl.expr <- as.matrix(read.table("dlbcl_expressiondata.csv",sep=","))
dlbcl.surv <- as.matrix(read.table("dlbcl_survdata.csv",sep=","))
dim(dlbcl.expr)
dim(dlbcl.surv)
dlbcl.surv
# Time, Censor, Genes
bigData <- cbind(dlbcl.surv,t(dlbcl.expr))


# divide data into training/test sets
set.seed(538)
train.index <- sample(c(1:240),120)
dlbcl.train.data <- dlbcl.expr[,train.index]
dlbcl.train.time <- dlbcl.surv[train.index,1]
dlbcl.train.cens <- dlbcl.surv[train.index,2]

dlbcl.test.data <- dlbcl.expr[,-train.index]
dlbcl.test.time <- dlbcl.surv[-train.index,1]
dlbcl.test.cens <- dlbcl.surv[-train.index,2]

dlbcl.expr.pvalue <- dlbcl.expr.zscore <- NA
for(i in 1:nrow(dlbcl.train.data)) {
  dlbcl.train.summary <- summary(coxph(Surv(dlbcl.train.time,dlbcl.train.cens) ~ dlbcl.train.data[i,]))
  dlbcl.expr.pvalue[i] <-dlbcl.train.summary$coef[5]
  dlbcl.expr.zscore[i] <- dlbcl.train.summary$coef[4]
}

min(dlbcl.expr.pvalue)
max(dlbcl.expr.pvalue)
hist(dlbcl.expr.pvalue,breaks=50)

# problem: expected number of false positives depends on number of features

# expected number of false positives (among truly null features): 
length(dlbcl.expr.pvalue)*.05

# expected number of false positives (among truly null features): 
length(dlbcl.expr.pvalue)*.01


# the false discovery rate assesses the proportion of false positives incurred
# when that particular test is called significant.
p.adjust(dlbcl.expr.pvalue,method="fdr")
dlbcl.expr.fdr <- p.adjust(dlbcl.expr.pvalue,method="fdr")

min(dlbcl.expr.fdr)
max(dlbcl.expr.fdr)
hist(dlbcl.expr.fdr,breaks=50)

ind <- sort.list(dlbcl.expr.pvalue)
dlbcl.expr.pvalue[ind]
write.table(cbind(dlbcl.expr.pvalue[ind],dlbcl.expr.fdr[ind]),
            "pvalue_fdr.csv",sep=",",quote=F,row.names=F,col.names=F)

sum(dlbcl.expr.pvalue < .15)
sum(dlbcl.expr.fdr < .15)

# current example: if pick p-value cutoff = 0.15
# implies 15% chance of false positive
# 7399 markers, expect about 1110 false positives if all features truly null
# in this data set, we have 1814 markers with p-value <= 0.15
# if all features truly null, expect 1110 of these to be false positives (61% of all identified markers)
# can't say what percent of these features are null without that assumption

# if pick fdr cutoff = 0.15
# expect 15% of all markers with fdr less than this to be false positives
# 13 with fdr value <= 0.15: expect 1-2 false positives among this 13



# fdr and p-values vs. z-scores
plot(as.numeric(dlbcl.expr.zscore),as.numeric(dlbcl.expr.pvalue),
     type="p",pch=19,cex=.3,col="blue",ylab="p-value (blue), fdr (red)",
     xlab="z-score",main="q-values and fdr vs. Z-scores")
lines(as.numeric(dlbcl.expr.zscore),as.numeric(dlbcl.expr.fdr),
      type="p",pch=19,cex=.3,col="red")


# fdr vs. p-values
plot(as.numeric(dlbcl.expr.pvalue),as.numeric(dlbcl.expr.fdr),
     type="p",pch=19,cex=.3,ylab="fdr",xlab="p-value",main="fdr vs. p-value")


# number of significant genes vs. fdr
numgenes <- NA
for(j in 1:length(unique(sort(dlbcl.expr.fdr)))) {
  numgenes[j] <- sum(sort(dlbcl.expr.fdr) <= unique(sort(dlbcl.expr.fdr))[j])
}
plot(unique(sort(dlbcl.expr.fdr)),numgenes,type="l",pch=19,cex=.3,ylab="number of significant genes",xlab="fdr",main="no. of significant genes vs. fdr")


# expected number of false positive genes vs. the total number of significant
#    genes given by the fdr
# i.e., number of markers declared signficant at a given fdr threshold
#    times that FDR threshold
plot(numgenes,unique(sort(dlbcl.expr.fdr))*numgenes,type="l",pch=19,cex=.3,ylab="number of expected false positives",xlab="number of significant genes",ylim=c(0,700),xlim=c(0,2000))





# plot of marker-by-marker p-values and fdr
plot(1:length(dlbcl.expr.fdr),as.numeric(dlbcl.expr.pvalue),type="p",pch=19,cex=.3,ylab="p-value (blue), q-value (red)",xlab="Marker",xaxt = "n",main="Expression Profiles Associated With Overall Survival",col="blue")
lines(1:length(dlbcl.expr.fdr),as.numeric(dlbcl.expr.fdr),type="p",pch=19,col="red",cex=.3)
abline(h=.15,lwd=2)

