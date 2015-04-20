# Build model using F-scale and past data
rm(list=ls())

# My Functions:
  rmCol  <- function(M,cn) M[,-c(which(is.element(colnames(M),cn)))]
  getCol <- function(M,cn) M[, c(which(is.element(colnames(M),cn)))]
  
# Clean Data:
  # Check for collinearity:
  # Remove repeating observations
  Dat <- read.csv("../Data/tornados.csv")
  Dat <-rmCol(Dat,"Year")
  Dat$Fscale <- as.factor(Dat$Fscale)
  # dim(dat) # 957 x 17

  #"Number"     "Month"      "Day"        "Date"       "Time"      
  #"State"      "Fscale"     "Injuries"   "Fatalities" "Loss"      
  #"CropLoss"   "StartLat"   "StartLon"   "EndLat"     "EndLon"    
  #"Length"     "Width"

  #dat <- subset(Dat,!duplicated(Dat$Number))
  dat <-Dat[-c(120,122,113,116,158,168,148,149,400,187,188,200,202,209,179,180,936),]


####################################################################
  # Random Forest:
  library(randomForest)
  set.seed(1)
  pos <- regexpr(":",dat$Time)
  dat$Time <- as.numeric(substr(dat$Time,1,pos-1))
  trainI <- sample(1:nrow(dat),800)

  library(foreach)
  library(doMC)
  registerDoMC(16)

  #small.dat <- rmCol(dat,c("Number","State","Date"))
  #small.dat <- rmCol(dat,c("Number","State","Date","Month","Day","Time",
  #                         "StartLat","StartLon","EndLat","EndLon"))
  #small.dat <- getCol(dat,c("Fscale","Injuries","Fatalities","Loss",
  #                          "CropLoss","Length","Width"))
  small.dat <- getCol(dat,c("Fscale","Month","Day","Time","Injuries",
                            "Fatalities","Loss","CropLoss","Length","Width"))

  one.it <- function(m) {
    forest <- randomForest(Fscale ~ .,data=small.dat,
                           subset=trainI,mtry=m,importance=T,ntree=1000)
    y.hat <- predict(forest,newdata=small.dat[-trainI,])
    err <- mean(y.hat!=dat$Fscale[-trainI])
    
    list("forest"=forest,"err"=err,"m"=m,"true"=dat$Fscale[-trainI],"pred"=y.hat)
  }
   
  result <- foreach(i=1:(ncol(small.dat)-2)) %dopar% one.it(i)
  forests <- lapply(result,function(x) x$forest)
  err.rate <- sapply(result,function(x) x$err)
  best.m <- which.min(err.rate)
  
  #forest <- randomForest(Fscale ~ .,data=small.dat,mtry=best.m,importance=T)
  forest <- randomForest(Fscale ~ .,data=small.dat,mtry=3,importance=T)
  # mtry=3=sqrt(P)
  importance(forests[[1]])
  varImpPlot(forests[[1]],main="Variable Importance")
  plot(forests[[1]])
  
  library(xtable)
  sink("../latex/out/confusion.tex")
    xtable(forest$confusion)
  sink()
  sink("../latex/out/output.txt")
    cat("\n"); paste("Error Rate:"); err.rate
    cat("\n"); paste("Best m:"); best.m
    cat("\n"); paste("Variable Importance:");importance(forests[[1]])
  sink()

 #TREE: #########################################################################
 library(tree)
 one.tree <- tree(Fscale ~ ., data=small.dat[trainI,])

 pdf("../latex/out/oneTree.pdf",width=20,height=10)
  plot(one.tree); text(one.tree)
 dev.off()
 pdf("../latex/out/forest.pdf")
  plot(forests[[1]])
 dev.off()
 pdf("../latex/out/importance.pdf")
  varImpPlot(forests[[1]],main="Variable Importance")
 dev.off()
