rm(list=ls())
options("width"=180)
gdp <- read.csv("GDP_data.csv")
library(glmnet)
grid <- 10^seq(-3,3,length=1000)
X <- model.matrix(GR6096 ~ .,data=gdp[,-c(1,2)])[,-1]
Y <- gdp$GR6096

this <- numeric(67)
this[1:67] <- FALSE

for( i in 1:67){
  if( length(unique(X[,i])) == 2) this[i] <- TRUE
} 

ind <- which(this==1)
x <- X[,]
x[,-ind] <- scale(X[,-ind])
x[,ind] <- X[,ind]
y <- scale(Y)

set.seed(1)
train <- sample(nrow(X),nrow(X)-10)
test <- (-train)
ridge.mod <- glmnet(x[train,],y[train],alpha=0,lambda=grid,standardize=F)

ridge_cv <- cv.glmnet(x[train,],y[train],nfolds=10,alpha=0,lambda=grid,standardize=F)
#551

pdf("ridge.pdf"); plot(ridge_cv); dev.off()

(best_lambda <- ridge_cv$lambda.min)

ridge_co <- predict(ridge.mod, s=best_lambda , type="coefficients")

ridge_pred <- predict(ridge.mod, s=best_lambda , newx=x[test,])
(mse <- mean((ridge_pred - y[test])^2))

x <- cbind(rep(1,60),x)

ridge_bet <-  solve(t(x) %*% x + diag(best_lambda,68)) %*% t(x) %*% y
k <- sum( diag( x %*% solve(t(x) %*% x + diag(best_lambda,68) ) %*% t(x) ))
n <- nrow(x)
s2 <- t(y - x %*% ridge_bet) %*% (y - x %*% ridge_bet) / (n-k)
std.err <- solve( t(x) %*% x + diag(best_lambda,68) ) %*% (t(x) %*% x) %*% solve( t(x) %*% x + diag(best_lambda,68) ) 
std.err <- sqrt(s2 * diag(std.err) )

ci_lo <- ridge_bet - qt(.975,n-k)* std.err
ci_hi <- ridge_bet + qt(.975,n-k)* std.err
ci <- cbind(ci_lo,ci_hi)

significant <- which(apply(ci,1,min)>0 | apply(ci,1,max) <0)

sort(abs(ridge_bet))
mse

#Tables:
t.val <- ridge_bet / std.err
p.val <- 2*pt(abs(t.val),nrow(x)-k,lower.tail=F) 

summary.model <- cbind(ridge_bet,ci_lo,ci_hi,std.err,t.val,p.val)
colnames(summary.model) <- c("Estimate","CI.Lo","CI.Hi","Std.Error","t-value","Pr(>|t|)") 
rownames(summary.model)[1] <- "(Intercept)"

summary.result <- summary.model[which(summary.model[,"Pr(>|t|)"] < .05),] #1
summary.result <- summary.result[order(summary.result[,"Pr(>|t|)"]),]
#summary.model[significant,] #2

library(xtable) 
xtab <- xtable(summary.result,digits=c(1,2,2,2,2,2,4))

sink("xtab.tex")
  xtab
sink()
#unlink("xtab.tex") #This removes the file "xtab.tex"
