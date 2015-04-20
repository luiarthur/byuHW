data<-read.table(header=T,text="
drugA drugB diff.VAS
0 0 -0.04
0 0 -0.02
0 0 -0.09
0 0 0.05
0 25 0.06
0 25 -0.18
0 25 0.00
0 25 -0.04
0 50 -0.36
0 50 -0.26
0 50 -0.27
100 0 0.00
100 0 -0.05
100 0 -0.26
100 0 0.10
100 25 -0.21
100 25 -0.06
100 25 -0.16
100 25 -0.10
100 50 -0.13
100 50 -0.03
100 50 -0.19
100 50 -0.30
150 0 -0.06
150 0 -0.17
150 0 -0.07
150 0 0.15
150 25 -0.26
150 25 0.05
150 25 -0.11
150 25 -0.09
")

data$drugA<-factor(data$drugA)
data$drugB<-factor(data$drugB)

model1<-lm(diff.VAS~drugA*drugB,data=data)
summary(model1)
anova(model1)

