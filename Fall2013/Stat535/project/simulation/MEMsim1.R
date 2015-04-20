# this simulation study investigates the effect of measurement error in X
# for the simple linear regression model

# True relationship:
#  mu_y = beta0 + beta1 * X
#
# But we observe
#  Y = mu_y + epsilon (as we have modeled)
#  W = X + u  (so we don't observe X --- we observe something close)

# number of simulated samples
nsim<-5000

# choose the values of the explanatory variable
#
#  Carroll, Ruppert, Stefanski Example 1.1.3 (Bioassay in a Herbicide Study)
#  suggests that while we intend to set each factor level in a designed experiment, 
#  the amount that actually is applied (or used in the experiment) is subject to variation.
#  (this is often referred to as the Berkson Model)
#
# consider a designed experiment with 5 replicates at each design point
#  (the x values represent dosages)
x<-rep(c(0,100,250,500),5)
n<-length(x)

# generate data from the probability model with parameters equal to
beta0 <- 5      # note: to approximate the power curve, we need to repeat simulation on a grid of beta0
                #       (for example, beta0=0, beta0=5, beta0=10, ... beta0=25)
beta1 <- 0.10   # note: to approximate the power curve, we need to repeat simulation on a grid of beta1
                #       (for example, beta1=0, beta1=0.01, beta1=0.02, ... , beta1=0.10)
sigma2 <- 100
# note: the levels of "measurement error effect" are
#          sigma2.w=0 (no effect)
#          sigma2.w=1000 (mild effect)
#          sigma2.w=10000 (strong effect)
sigma2.w <- 10000


# simulation
hat.beta0<-rep(0,nsim)
hat.beta1<-rep(0,nsim)
s2<-rep(0,nsim)
L0<-rep(0,nsim)
U0<-rep(0,nsim)
L1<-rep(0,nsim)
U1<-rep(0,nsim)
pval0<-rep(0,nsim)
pval1<-rep(0,nsim)
for(sim in 1:nsim){

    # observed data
    y<-beta0+beta1*x+sqrt(sigma2)*rnorm(n)
    w<-x+sqrt(sigma2.w)*rnorm(n)

    # compute regression estimates
    ybar<-mean(y)
    wbar<-mean(w)
    sum.w2<-sum(w^2)
    sum.wy<-sum(w*y)
    hat.beta1[sim]<-(sum.wy - n*wbar*ybar)/(sum.w2 - n*wbar^2)
    hat.beta0[sim]<-ybar-hat.beta1[sim]*wbar
    s2[sim]<-sum((y-hat.beta0[sim]-hat.beta1[sim]*w)^2)/(n-2)
    # compute CIs
    L0[sim]<-hat.beta0[sim] - qt(0.975,n-2)*sqrt(s2[sim]*((1/n)+(wbar^2/(sum.w2 - n*wbar^2))))
    U0[sim]<-hat.beta0[sim] + qt(0.975,n-2)*sqrt(s2[sim]*((1/n)+(wbar^2/(sum.w2 - n*wbar^2))))
    L1[sim]<-hat.beta1[sim] - qt(0.975,n-2)*sqrt(s2[sim]/(sum.w2 - n*wbar^2))
    U1[sim]<-hat.beta1[sim] + qt(0.975,n-2)*sqrt(s2[sim]/(sum.w2 - n*wbar^2))
    # compute pvalues for HT
    test.t<-hat.beta0[sim]/sqrt(s2[sim]*((1/n)+(wbar^2/(sum.w2 - n*wbar^2))))
    pval0[sim]<-2*(1-pt(abs(test.t),n-2))
    test.t<-hat.beta1[sim]/sqrt(s2[sim]/(sum.w2 - n*wbar^2))
    pval1[sim]<-2*(1-pt(abs(test.t),n-2))

}

# Results for the specified parameter values:

# Bias
print(c(mean(hat.beta0-beta0),mean(hat.beta1-beta1),mean(s2-sigma2)))

# MSE
print(c(mean((hat.beta0-beta0)^2),mean((hat.beta1-beta1)^2),mean((s2-sigma2)^2)))

# 95% Confidence Interval Coverage for beta0
print(mean(L0<beta0 & beta0<U0))
# 95% Confidence Interval Coverage for beta1
print(mean(L1<beta1 & beta1<U1))

# Power of Ho: beta0=0 (at the specified value of beta0)
print(mean(pval0<0.05))
# Power of Ho: beta1=0 (at the specified value of beta1)
print(mean(pval1<0.05))
