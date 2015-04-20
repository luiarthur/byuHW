# My own function :)
se.pcnt <- function(p,km,e=.05){# p is a percentage for a percentile
  # Computes the standard error of a given percentile
  se.S <- km$std.err[which(km$surv<p)[1]]
  perc <-  km$time[which(km$surv<p)[1]]
  u <- tail(km$time[which(km$surv >= 1-p+e)],1)
  l <- km$time[which(km$surv <= 1-p-e)][1]
  S.u <- tail(km$surv[which(km$surv >= 1-p+e)],1)
  S.l <- km$surv[which(km$surv <= 1-p-e)][1]

  f <- (S.u-S.l) / (l-u)
  se <- se.S / f

  CI <- matrix(c(perc,qnorm(c(.025,.975),perc,se)),1,3)
  colnames(CI) <- c("Estimate","CI.Lower","CI.Upper")
  list("se"=se,"percentile"=perc,"CI"=CI)
}


