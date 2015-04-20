rm(list=ls())
gdp <- read.csv("GDP_data.csv")
library(glmnet)
grid <- 10^(seq(-5,0,length=1000))
X <- model.matrix(GR6096 ~ .,data=gdp[,-c(1,2)])[,-1]
Y <- gdp$GR6096

set.seed(1)
train <- sample(nrow(X),nrow(X)-10)
test <- (-train)
ridge.mod <- glmnet(X[train,],Y[train],alpha=0,lambda=grid)
ridge_cv <- cv.glmnet(X[train,],Y[train],nfolds=10,alpha=0,lambda=grid)

plot(ridge_cv)

(best_lambda <- ridge_cv$lambda.min)

ridge_co <- predict(ridge.mod, s=best_lambda , type="coefficients")
ridge_pred <- predict(ridge.mod, s=best_lambda , newx=X[test,])

(mse <- mean((ridge_pred - gdp$GR6096[test])^2))

this <- numeric(67)
this[1:67] <- FALSE

for( i in 1:67){
  if( length(unique(X[,i])) == 2) this[i] <- TRUE
} 

ind <- which(this==1) # ind = index for categorical variables
x <- X[,]
x[,-ind] <- scale(X[,-ind])
x[,ind] <- X[,ind]
x <- cbind(1,x)
y <- scale(Y)
ridge_bet <-  solve(t(x) %*% x + diag(best_lambda,68)) %*% t(x) %*% y
q <- solve ( t(x) %*% x + best_lambda * diag(ncol(x)) )
k <- sum(diag(x %*% q %*% t(x)))

s2 <- t(y - x %*% ridge_bet) %*% (y - x %*% ridge_bet) / (nrow(x)-k)
s2 <- as.numeric(s2)
sort(abs(ridge_bet))


se <- sqrt(diag(s2*q %*% t(x)%*%x * q))
t <- ridge_bet / se
p.value <- 2*pt(abs(t),nrow(x)-k,lower.tail=F)

summary.model <- cbind(ridge_bet,se,t,p.value)
colnames(summary.model) <- c("Estimate","Std. Error","t-value","Pr(>|t|)") 
rownames(summary.model)[1] <- "(Intercept)"

head(summary.model)
# Results, Conclusions
# SIGNIFICANT:
# Intercept, Catholic, CIV72, CIVLIBTERTY, CONFUCIUS, DENS, EAST, EUROPE,
# 
