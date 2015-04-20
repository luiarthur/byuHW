health <- rep(c("support","oppose"),4)
info <- rep(rep(c("support","oppose"),each=2),2)
gender <- c(rep("M",4),rep("F",4))
y <- c(76,160,6,25,114,181,11,48)

GH.GI <- glm(y ~ gender+info+health+gender:health+
                 gender:info,family=poisson(link="log"))
GH.HI <- glm(y ~ gender+info+health+gender:health+
                 health:info,family=poisson(link="log"))
GI.HI <- glm(y ~ gender+info+health+gender:info+
                 health:info,family=poisson(link="log"))
GH.GI.HI <- glm(y ~ gender*info*health-gender:info:health,
                    family=poisson(link="log"))

#a)
GH.GI$deviance #HIGHEST
GH.HI$deviance
GI.HI$deviance
GH.GI.HI$deviance

#b)
summary(GH.GI.HI)
gh <- -.2516 + c(-1,1)*1.96*.1749
exp(gh)

gi <- .4636 + c(-1,1)*1.96*.2406
exp(gi)
