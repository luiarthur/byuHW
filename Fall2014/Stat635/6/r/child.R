library(lme4)

dat <- read.table("../child.dat")
kid <- dat$V1
sex <- dat$V2
age <- dat$V3
len <- dat$V4

### Mixed effects model, accounting for random effects for each individual
mod_me = lmer(len ~ age+sex+age*sex+(1|kid))


### Linear model, not accounting for random effects for each individual
mod_0 = lm(len ~ age+sex+age*sex)

s_me = summary(mod_me)
s_0 = summary(mod_0)

plot(s_0$resid)
points(s_me$resid, pch=20, col='red')

### Take aways:
# y = Xb + Zu + e.
# b is fixed
# u ~ Normal(0, v*I)
# e ~ Normal(0, sig2* I)
# We use random effects model when we have repeated measures on 
# individuals that are assumed to be independent. In this case, we 
# model kids elbow lengths, with regressors: sex and age. We have repeated
# measures for kids. The measurements for a kid are surely correlated,
# but kid1 and kid2 are independent. We can reasonably use a random intercept
# for each kid. This reduces the residuals of the model, and the coefficients
# of the fixed effects will have smaller variances. The point estimates 
# of the fixed effects will remain the same in the model that doesn't
# account for different kids. For predictions for a new kid, just use the fixed
# effects. For predictions on one of the experimented kids, use the random
# effects parameters also. Finally, the parameter u is random; not fixed.
# This ensures the parameters have mean 0. This is a constraint. In this random
# intercept model, without a prior mean, the model is unidentifiable. This is 
# because we also have an overall intercept in the fixed effects.  Without the
# prior mean of 0, we have a non-full rank design matrix. 
# If we do everything correctly, we'll have a more natural interpretation of the 
# problem.
