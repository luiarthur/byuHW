rm(list=ls())
options("width"=180)
GDP <- read.csv("GDP_data.csv",sep=",")

get.estimate <- function(gdp=GDP,z=1){
  # Exploring the data:  

    #any(is.na(gdp))
    
    #pairs(gdp[,1:10])
    #pairs(gdp[,c(3,11:20)])
    #pairs(gdp[,c(3,21:30)])
    #pairs(gdp[,c(3,31:40)])
    #pairs(gdp[,c(3,41:50)])
    #pairs(gdp[,c(3,51:60)])
    #pairs(gdp[,c(3,61:70)])

    #vis <- gdp[,c("CODE","GR6096")]
    #vis <- vis[order(vis[,2]),]
    #plot(vis,las=2,cex.axis=.6)

    #GR <- vis[,2]
    #plot(GR,pch=20)

  # For TESTING:
  # gdp <- GDP

  X <- model.matrix(GR6096 ~ .,data=gdp[,-c(1:2)])[,-1]
  Y <- gdp$GR6096

  library(glmnet)
  grid <- 10^seq(10,-2,length=100)
  ridge.mod <- glmnet(X,Y,alpha=0,lambda=grid) # alpha=0 => rigde,
                                               # alpha=1 => lasso.
  # Exploring:
  #  plot(ridge.mod)
  #  coef(ridge.mod)[,2]
  #  dim(coef(ridge.mod))
  #  ridge.mod$lambda

  #temp <- predict(ridge.mod,s=50, type="coefficients")  

  #set.seed(1)

  test.size <- 10 
  train <- sample(1:nrow(X),nrow(X)-test.size)
  test <- (-train)
  cv.out <- cv.glmnet(X[train,], Y[train], alpha=0, nfolds=10) #ridge
  #plot(cv.out)
  bestlam <- cv.out$lambda.min

  ridge.pred <- predict(ridge.mod, s=bestlam, newx=X[test,])
  mse <- mean((ridge.pred - Y[test])^2) #Print
  #cbind(ridge.pred, Y[test])

  out <- glmnet(X,Y,alpha=0)
  model <- as.matrix(predict(out, type="coefficients", s=bestlam))
  #sorted.model <- as.matrix(model[order(abs(model[,1]),decreasing=T),])
  colnames(model) <- "Estimate"
  
  if (z%%20==0) cat(">")

  return(model)
}

# Estimate, S.E., t-value Pr(>|t|)
# top.factors <- head(sorted.model,10)

library(foreach)
library(doMC)
registerDoMC(16)

get.se <- function(B=1000){
  cat(paste(rep("#",50),collapse="")); cat("\n")
  K <- nrow(GDP) # How do I choose K?
  boot <- GDP[sample(1:nrow(GDP), K, replace=T),]
  estimates <- foreach(b=1:B,.combine=cbind) %dopar% get.estimate(boot,b)
  beta.se <- apply(estimates,1,sd)
  cat("\n")
  return(beta.se)
}

Estimates <- get.estimate(GDP)
SE <- get.se(1000)
t <- Estimates/SE
p.value <- NA # 2*pt(t,length(Estimates)-1)

summary.model <- cbind(Estimates,SE,t,p.value)
colnames(summary.model) <- c("Estimate","Std. Error","t-value","Pr(>|t|)")

sorted.model <- as.matrix(summary.model[order(abs(summary.model[,1]),
                          decreasing=T),])
head(sorted.model)
