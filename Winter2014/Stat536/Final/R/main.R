# Code to be run in UNIX/linux environment
system("mkdir -p ../latex/raw")
rm(list=ls())
library(xtable)

# Parameter Declarations:
d <- 3 # df for natural spline

# My Functions:
  rmCol  <- function(M,cn) M[,-c(which(is.element(colnames(M),cn)))]
  getCol <- function(M,cn) M[, c(which(is.element(colnames(M),cn)))]
  
  case <- function(test,a,b) {
    if (test==T) {
      a
    } else {
      b
    }
  }

# Read & Clean Data:
  Dat <- read.csv("../Data/Tulips.csv")
  dat <- rmCol(Dat,c("YearCollected","DayCollected"))
  dat <- dat[-which(dat$Population==12),]
  dat <- dat[order(dat$Pop,dat$Chill),]
  #dat$Germinated <- ifelse(dat$Germinated=="Y",1,0)
  dat$Population <- as.factor(dat$Population)
  colnames(dat) <- c("Pop","Chill","Germ")
  counts <- table(dat)
  props  <- counts[,,2] / (counts[,,1] + counts[,,2])
  #colnames(props) <- paste("Pop",0:6 * 2)
  #rownames(props) <- paste(1:11,"Weeks")
  round(props,2)
  
  #sink("../latex/raw/props.tex")
  #  xtable(props)
  #sink()  

# Exploratory Plots  
  pop <- lapply(as.list(1:11), function(x) dat[which(dat$Pop==x),])

  #pdf("../latex/raw/pairs.pdf")
  #  pairs(dat,col="#9999ff22",pch=20,cex=3,main="Overall")
  #dev.off()
  #for (i in 1:11) {
  #  pairs(pop[[i]],col="#9999ff22",pch=20,cex=3,main=paste("Pop",i))
  #}  

  plot.dat <- function(i=0,add=F) {
    do <- case(add,lines,plot)
    d  <- case(i==0,apply(props,2,mean),props[i,])
    a  <- case(i==0,"All Populations",paste("Popultation",i))
    x  <- seq(0,12,by=2)
    do(x,d,type='p',ylim=c(0,1),main=paste("Germination Rates for",a),pch=20,
       ylab="Probability",xlab="Chill Time (Weeks)")
  }

  pdf("../latex/raw/dat.pdf")
    par(mfrow=c(6,2),mar=rep(3,4))
    for(i in 0:11) plot.dat(i)
    par(mfrow=(c(1,1)))
  dev.off()

# Natural Spline:
  library(splines)
  #create Full Model:
    temp <- dat
    temp$Pop <-as.factor(0)
    bigD <- rbind(dat,temp)
    
    cv.fullMod <- function(i) {
      ind <- ((i-1)*210+1):case(i==12,nrow(bigD),i*210)
      bd <- bigD[ind,]
      n <- length(ind)
      trainI <- sample(1:n,round(n*.8))
      fm <- glm(Germ~ns(Chill,d),data=bd[trainI,],family=binomial)
      pred <- predict(fm,list("Chill"=bd$Chill[-trainI]),type="response") 
      G <- ifelse(pred<.5,"N","Y")
      mean(bd$Germ[-trainI] != G)
    }
    
    B <- 100
    get.err.rate <- function(y) apply(matrix(1:12),1,function(x) cv.fullMod(x))
    err.rate <- lapply(as.list(1:B),get.err.rate)
    err.r <- t(as.matrix(Reduce("+",err.rate) / B))
    colnames(err.r) <- c(paste("Population",1:11),"Overall")

    #sink("../latex/raw/err.tex")
    #  xtable(err.r)
    #sink()

    mod.all <- glm(Germ~ns(Chill,d),data=dat,family=binomial)
    fullMod <- glm(Germ~ns(Chill,d)+Pop+ns(Chill,d)*Pop,data=bigD,family=binomial)
    mod <- lapply(as.list(1:11),function(x) glm(Germ~ns(Chill,d),data=pop[[x]],
                                                family=binomial))
    plot.pred <- function(i=0,compare=T) {
      x0 <- seq(0,12,length=1000)
      pred <- case(i==0,
                    predict(mod.all,list("Chill"=x0),type="response"),
                    predict(mod[[i]],list("Chill"=x0),type="response")
                  )
      a <- case(i==0,"All Populations",paste("Popultation",i))
      plot(x0,pred,type='l',ylim=c(0,1),main=paste("Germination Rates for",a),
           col="purple",lwd=3)
      if (compare) plot.dat(i,add=T)
    }
    
    plot.all <- function() {
      x0 <- seq(0,12,length=1000)
      par(mfrow=c(6,2),mar=rep(3,4))
      for (i in 0:11) {
        plot.pred(i) 
        #temp <- predict(fullMod,newdata=list("Chill"=x0,
        #        "Pop"=as.factor(rep(i,1000))),type="response")
        #lines(x0,temp,ylim=c(0,1),type='l',col='red')
      }  
      par(mfrow=c(1,1))
    }
    
    pdf("../latex/raw/all.pdf")
      plot.all()
    dev.off()

  one.compare <- function(i,j) {
    pops <- rbind(pop[[i]],pop[[j]])
    fm <- glm(Germ ~ Pop+ns(Chill,d)+Pop*ns(Chill,d),data=pops,family=binomial)
    rm <- glm(Germ ~ ns(Chill,d),data=pops,family=binomial)
    pval <- anova(rm,fm,test="Chisq")$Pr[2]
    pval
  }
  
  compare.all <- function() {
    comp <- matrix(0,11,11) 
    for (i in 1:11) {
      for (j in setdiff(1:i,i)) {
        comp[i,j] <- one.compare(i,j)
      }
    }
    ind <- which( (comp > .05 / (choose(11,2))) & (comp!=0) ) # Bonferroni
    pairs <- paste("(",ifelse(ind%%11==0,11,ind%%11),",",
                   ifelse(ind%%11==0,ind%/% 11,ind%/%11+1),")",
                   sep="",collapse=", ")
    list("M"=comp,"pairs"=pairs)
  }

  options("width"=180)
  comp <- compare.all()
  M <- comp$M
  same <- comp$pairs
  options("width"=80)

  sink("../latex/raw/comparisons.txt")
    M
    same
  sink()  
  # 1:
  # The effect of chilling time is not the same across different populations.
  # The following populations behave the most similiarly under different chill times: # (3,2), (4,2), (10,4), (7,6), (10,6), (11,6), (10,7), (11,7), (11,10)
  
  boot <- function(fn,i,B=10000) {
    N <- nrow(dat)
    n <- nrow(pop[[1]])

    get.mod <- function() {
      trI.m.all <- sample(1:N,N,replace=T)
      trI.m <- lapply(as.list(1:11), function(x) sample(1:n,n,replace=T) )

      p <- lapply(as.list(1:11), function(x) dat[which(dat$Pop==x),][trI.m[[x]],])

      m.all <- glm(Germ~ns(Chill,d),data=dat[trI.m.all,],family=binomial)
      m <- lapply(as.list(1:11),function(x) glm(Germ~ns(Chill,d),data=p[[x]],
                                                  family=binomial))
      list(m,m.all)
    }

    doit <- function(i) {
      mds <- get.mod()
      m.all <- mds[[2]]
      m <- mds[[1]]
      
      fn(i,m,m.all)
    }

    library(doMC)
    registerDoMC(16)

    btstrp <- foreach(b=1:B,.combine=rbind) %dopar% doit(i)
    se <- sd(btstrp)
    est <- mean(btstrp)
    out <- matrix(c(est,qnorm(c(.025,.975),est,se)),1,3)
    out
  } 

  # 2:
  # Ideal Chilling Time:
  best.chill.time <- function(i=12,md=mod,md.all=mod.all) {
    x0 <- seq(0,12,length=1000)
    pred <- case(i==12,
                  predict(md.all,list("Chill"=x0),type="response"),
                  predict(md[[i]],list("Chill"=x0),type="response")
                )
    temp <- pred[which.max(pred)]
    name <- as.numeric(names(temp))
    best.chill.time <- x0[name]
    best.chill.time
  }


  # Answer 2: Do I need Uncertainties?
  # Best across populations is 12 
  # Best varies by population
  best.chill.times <- t(apply(matrix(1:12),1,function(x) t(boot(best.chill.time,x,100))))
  colnames(best.chill.times) <- c("Estimate","95% Lower.CI","95% Upper.CI")
  rownames(best.chill.times) <- c(paste("Population",1:11),"All Populations")

  library(xtable)
  sink("../latex/raw/best.chill.tex")
    xtable(best.chill.times)
  sink()  

  plot.order <- function(m,ylb,xlb,mn,ylm,eps,yax) {
    y <- m[-12,]
    i <- rev(order(y[,1]))
    y <- y[i,]
    x <- 1:11
    plot(y[,1],pch=20,col="blue",xlim=c(1,11),ylim=ylm,
         ylab=ylb,xlab=xlb,
         main=mn,axes=F)
    axis(side=1,at=1:11,i) 
    axis(side=2,at=yax) 
    segments(x,y[,2],x,y[,3])
    epsilon <- eps
    segments(x-epsilon,y[,2],x+epsilon,y[,2])
    segments(x-epsilon,y[,3],x+epsilon,y[,3])
  }

  pdf("../latex/raw/bct.pdf")
    plot.order(best.chill.times,xlb="Population",ylb="Best Chilling Time (Weeks)",
               mn="Best Chilling Times",ylm=c(0,16),eps=.2,yax=0:16)
  dev.off()  
  ###############################################################################

  # 3: What effect will a decrease from 10 to 8 weeks of
  #    chilling time have for tulips?
  effect <- function(i,md=mod,md.all=mod.all) {#Effect of changing time from 10 to 8
    m <- case(i==12,md.all,md[[i]])
    x0 <- c(8,10);  x0 <- as.data.frame(x0); colnames(x0) <- "Chill"
    pred.p <- predict(m,x0,type="response")
    pred.p[1] - pred.p[2]
  }
  
  # Answer 3: Do I need uncertainties?
  # decrease by  -0.04123726  globally
  # the change varies
  effect.10.to.8 <- t(apply(matrix(1:12),1,function(x) t(boot(effect,x,100))))

  colnames(effect.10.to.8) <- c("Estimate","95% CI.Lo","95% CI.Hi")
  rownames(effect.10.to.8) <- c(paste("Population",1:11),"All Populations")

  sink("../latex/raw/effect.10.to.8.tex")
    xtable(effect.10.to.8)
  sink()  

  pdf("../latex/raw/eff.pdf")
    plot.order(effect.10.to.8,xlb="Population",ylb="Effect (Weeks)",ylm=c(-.8,.2),
               mn="Effect of Decrease in Chilling Time from 10 to 8 Weeks",eps=.2,
               yax=seq(-.8,.2,by=.1))
    abline(h=0,col='red')
  dev.off()  

