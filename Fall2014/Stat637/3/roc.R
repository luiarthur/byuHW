#8
dat <- read.table("../2/case.txt",header=T)
y <- dat$Rem
m <- dat$Cases
x <- dat$L1

#a)
mod <- glm(y/m~x,weights=m,family=binomial(link="logit"))
pred <- predict(mod,type="response")

calc.sens.spec <- function(n.p0=1e3) {
  p0.seq <- seq(0,1,len=n.p0)

  O <- matrix(0,n.p0,2)
  colnames(O) <- c("1-Specificity","Sensitivity")

  for (i in 1:n.p0) {
    p0 <- p0.seq[i]
    p.i <- ifelse(y/m>p0,1,0)
    p.i.hat <- ifelse(pred>p0,1,0)

    true.ind <- which(p.i==1)
    false.ind <- which(p.i==0)
    sens <- mean(p.i[true.ind]==p.i.hat[true.ind])
    spec <- mean(p.i[false.ind]==p.i.hat[false.ind])

    O[i,1] <- 1-spec
    O[i,2] <- sens
  }

  O
}

M <- calc.sens.spec()
plot(M,type="o",col="gray30",lwd=2,main="ROC for Logit Model");abline(0,1)

# Comment on the model fit according to these diagnostics and discuss how this
# changes or validates your prior conclusions about the model.
# The model predicts better than guessing. This means that the logit model was
# a good model choice.
