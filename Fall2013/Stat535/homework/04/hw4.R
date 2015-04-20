n <- 50
x <- c(1,3,5,2,4); k <- length(x)
x <- rep(jitter(X), n) * (1:n)
X <- matrix(x,n,k)

Ys <- ( Y - mean(Y) ) / sd(Y)
Xs <- ( X - apply(X,2,mean) ) / apply(X,2,sd)

mod <- lm( Y ~ X )
modC <- lm( cY ~ cX )
