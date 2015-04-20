rm(list=ls())

# Data Readin and Cleaning:
  cars <- read.csv("Cars.csv",header=T)
  cars$cc[81] <- 1600 # Because this was obviously a mistake
  cars <- cars[,-c(1,2,4,13)] #1:  Id removed because it's an index
                              #2:  Model removed because all Corolla
                              #4:  Age removed because redundant with year
                              #13: Cylinders removed because all are 4 Cyls
                              
  cols <- (1:ncol(cars))[-c(1,3,10)] # Price, Miles, Weight are quantitative
  for (i in cols) {
    cars[,i] <- as.factor(cars[,i])  # Every other covariate is a factor
  }


# Spline for miles and other variables are categorical:
  mod <- smooth.spline(cars$Miles,cars$Price,cv=T)
  lambda <- mod$lambda

  plot.smooth.spline <- function(){
    pdf("priceWeight.pdf")
    plot(cars$Miles,cars$Price,pch=20,cex=.7,
         xlab="Price",ylab="Miles",
         main="Price Vs. Miles")

    lines(mod,lwd=3,col='blue')
    legend("topright",legend="Smoothing Spline",
           col="blue",lwd=3)
    dev.off()       
  }
  
  #plot.smooth.spline()

# GAM:
  #library(splines)
  #library(gam)
  #cars.gam <- cars; cars.gam$Miles <- s(cars$Miles)
  #gam.mod <- gam(Price ~ ., data=cars.gam)
  #summary(gam.mod)
  #plot.gam(gam.mod, se=T, col='blue', ask=T)

# GAM: mgcv
  #library(gam)
  library(mgcv)
  temp <- cars[,2]; ctemp <- colnames(cars)[2]
  cars[,2] <- cars[,3]; colnames(cars)[2] <- colnames(cars)[3]
  cars[,3] <- temp; colnames(cars)[3] <- ctemp

  #form <- paste(colnames(cars[,-c(1,2)]),collapse="+")
  #form <- paste("Price ~", "s(cars$Miles)+", form)
  #gam.mod <- gam(as.formula(form), data=cars)
  #sum.gam.mod <- summary(gam.mod)

  form <- paste(colnames(cars[,-c(1,2,6,9,12:14,17,18,20,21)]),collapse="+") # Because Colors, Doors,
                                                                 # ABS, AirBag, 
                                                                 # Boardcomputer, CD_Player,
                                                                 # Power_Steering, Radio insignificant
  form <- paste("Price ~", "s(cars$Miles)+", form)

  gam.mod <- gam(as.formula(form), data=cars)
  sum.gam.mod <- summary(gam.mod)

  #low.mod <- gam(Price~1,data=cars)
  #forw.gam.mod <- step(low.mod,scope=list(lower=low.mod,upper=gam.mod),upper=
  #                     gam.mod,direction="both")

  # Code to get P.I.:
  #pred   <- predict.gam(gam.mod, newdata=cbind(cars$Miles[testI],cars[testI,-c(1,3)]), se.fit=T)
  #X <- data.frame(model.matrix(Price ~ ., data=cars)[,-c(6,9,12:14,17,18,20,21)])
  #pred    <- predict.gam(gam.mod, newdata=cars[1:100,-1], se.fit=T, na.action=na.omit, family="gaussian")
  #testI  <- sample(1:nrow(cars),100,replace=T)

  #new.cars <- data.frame(cars[,-c(6,9,12:14,17,18,20,21)])
  #X <- data.frame(model.matrix(Price ~ ., data=new.cars))[,-1]
  #pred    <- predict.gam(gam.mod, new.cars[,-1], se.fit=T, na.action=na.omit)
  pred    <- predict.gam(gam.mod, se.fit=T, na.action=na.omit)
  pred.se <- sqrt(pred$se.fit^2+gam.mod$sig2)
  pi.low  <- pred$fit - qt(.975,df=gam.mod$df.residual) * pred.se
  pi.up   <- pred$fit + qt(.975,df=gam.mod$df.residual) * pred.se

  coverage  <- mean(pi.low < cars$Price & cars$Price < pi.up)
  mean.pred <- mean(pred$fit)
  mean.PI   <- mean(pi.up-pi.low)

  plot.resid <- function(){
    pdf("resid.pdf")
    plot(gam.mod$residuals)
    dev.off()
  }  
  plot.qqnorm <- function(){
    pdf("qq.pdf")
    qqnorm(gam.mod$residuals)
    dev.off()
  }

  library(xtable)
  xtab.p <- xtable(sum.gam.mod$p.table)
  sink("xtabP.tex"); xtab.p; sink()
  xtab.s <- xtable(sum.gam.mod$s.table)
  sink("xtabS.tex"); xtab.s; sink()
  
  plot.resid(); plot.qqnorm(); plot.smooth.spline()
