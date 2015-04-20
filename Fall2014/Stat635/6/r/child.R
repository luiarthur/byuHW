library(lme4)

dat <- read.table("../child.dat")
kid <- dat$V1
sex <- dat$V2
age <- dat$V3
len <- dat$V4

lmer(len ~ age+sex+age*sex+(1|kid))
