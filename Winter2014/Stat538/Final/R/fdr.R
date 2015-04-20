# FDR: False Discovery Rate:= false positive / (false positive + true positive)

# zenburn scheme for linux?
# Questions:
#   1) How do I know if 1 means death or censor?
#   2) Do I use overall survival?

rm(list=ls())
# Load Libraries:
  library(rms)
  source("../../MySurv/mysurv.R",chdir=T)

# Clean Data:
  # Expression Data:
  exprss <- read.csv("../Data/BladderCancer_expression.csv")
  # Make the first column the column name
  rownames(exprss) <- exprss[,1]
  exprss <- exprss[,-1]

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

  sink("out/results.txt")
    cens.rate <- mean(clinic$cens==0)
    cat(paste("Number of Obs: ",N)); cat("\n")
    cat(paste("Number of Obs Censored: ",sum(clinic$cens==0))); cat("\n")
    cat(paste("Censoring Rate:",cens.rate)); cat("\n")
    cat(paste("Number of Genes: ",P)); cat("\n")
    cat("\n")
  sink() 

# Exploration:
  library(rms)
  KM <- survfit(Surv(time,cens) ~ 1, type="kaplan-meier",data=as.data.frame(bigD))
  #pdf("latex/raw/KM.pdf")
  #  survplot(KM)
  #dev.off()  
  #se.pcnt(.5,KM)
  #bigD[1:5,1:5] # row are subjects, columns are variables

# FDR:
  # inspect all code below here: 31 March, 2014.
  # Divide Data into Training/Test sets
    set.seed(538)
    train.ind <- sample(c(1:N),round(N/2))
    train.set <- bigD[train.ind,]
    test.set  <- bigD[-train.ind,]

    library(foreach)
    library(doMC)
    registerDoMC(16)

    one.summary <- function(i) {
      train.summary <- summary(coxph(Surv(train.set$time,train.set$cens) ~ train.set[,i]))
      p <- train.summary$coef[5]
      z <- train.summary$coef[4]
      matrix(c("pval"=p,"zscore"=z),1,2)
    }

    system.time(
      pz <- foreach(i=1:(ncol(train.set)-2),.combine=rbind) %dopar% one.summary(i) # Takes 35s
    )  
    colnames(pz) <- c("pval","zscore")
    pz <- data.frame(pz)

    ############################
    # Need to sink results!!!!!#
    ############################

    min(pz$pval)
    max(pz$pval)
    hist(pz$pval,breaks=50) # Need to see this. Haven't seen this yet.
    # The histogram indicates that out Genes are NOT SIGNIFICANTLY DIFFERENT?

    ###############################################################################
    # problem: expected number of false positives depends on number of features

    # expected number of false positives (among truly null features): 
    length(pz$pval)*.05

    # expected number of false positives (among truly null features): 
    length(pz$pval)*.01


    # the false discovery rate assesses the proportion of false positives incurred
    # when that particular test is called significant.
    fdr <- p.adjust(pz$pval,method="fdr")

    min(fdr)
    max(fdr)
    hist(fdr,breaks=50)

    ind <- sort.list(pz$pval)
    #pz$pval[ind]
    write.table(cbind(pz$pval[ind],fdr[ind]),"out/pval_fdr.csv",sep=",",quote=F,row.names=F,col.names=F)

    # Choose threshold. Maybe not .15?
    thresh <- .15        #.15 = original
    sum(pz$pval < thresh)
    sum(fdr     < thresh)
    ##########################################################################

    # fdr and p-values vs. z-scores
    plot(pz$zscore,pz$pval,
         type="p",pch=19,cex=.3,col="blue",ylab="p-value (blue), fdr (red)",
         xlab="z-score",main="q-values and fdr vs. Z-scores")
    lines(pz$zscore,fdr,type="p",pch=19,cex=.3,col="red")

    ## fdr vs. p-values
    plot(as.numeric(pz$pval),as.numeric(fdr),
         type="p",pch=19,cex=.3,ylab="fdr",xlab="p-value",main="fdr vs. p-value")


    # number of significant genes vs. fdr
    numgenes <- NA
    for(j in 1:length(unique(sort(fdr)))) {
      numgenes[j] <- sum(sort(fdr) <= unique(sort(fdr))[j])
    }
    plot(unique(sort(fdr)),numgenes,type="l",pch=19,cex=.3,ylab="number of significant genes",xlab="fdr",main="no. of significant genes vs. fdr")


    # expected number of false positive genes vs. the total number of significant
    #    genes given by the fdr
    # i.e., number of markers declared signficant at a given fdr threshold
    #    times that FDR threshold
    plot(numgenes,unique(sort(fdr))*numgenes,type="l",pch=19,cex=.3,ylab="number of expected false positives",xlab="number of significant genes",ylim=c(0,700),xlim=c(0,2000)) # Why 2000?



    # plot of marker-by-marker p-values and fdr
    # fdr is q-value?
    plot(1:length(fdr),pz$pval,type="p",pch=19,cex=.3,ylab="p-value (blue), q-value (red)",xlab="Marker",xaxt = "n",main="Expression Profiles Associated With Overall Survival",col="blue")
    lines(1:length(fdr),fdr,type="p",pch=19,col="red",cex=.3)
    abline(h=thresh,lwd=2)

