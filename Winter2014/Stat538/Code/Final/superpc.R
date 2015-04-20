library(survival)
library(superpc)


##################################################
### read in data

dlbcl.data <- as.matrix(read.table("dlbcl_expressiondata.csv",sep=","))
dlbcl.surv <- as.matrix(read.table("dlbcl_survdata.csv",sep=","))
dim(dlbcl.data)
dlbcl.surv
dlbcl.time <- dlbcl.surv[,1]
dlbcl.cens <- dlbcl.surv[,2]

# divide data into training/test sets
set.seed(3)
train.index <- sample(c(1:240),120)
dlbcl.train.data <- dlbcl.data[,train.index]
dlbcl.train.time <- dlbcl.surv[train.index,1]
dlbcl.train.cens <- dlbcl.surv[train.index,2]

dlbcl.test.data <- dlbcl.data[,-train.index]
dlbcl.test.time <- dlbcl.surv[-train.index,1]
dlbcl.test.cens <- dlbcl.surv[-train.index,2]

nrow(dlbcl.train.data) # markers
ncol(dlbcl.train.data) # subjects

### can't do PCA when p > n: prcomp just ignores any p that are beyond n
### stdev matrix should be a pxp matrix
res.prcomp <- prcomp(t(dlbcl.train.data),center=T,scale=T,retx=T)
res.prcomp$sdev

res.prcomp <- prcomp(dlbcl.train.data,center=T,scale=T,retx=T)
res.prcomp$sdev

###############################################
#### superpc ###############

featurenames <- paste("feature",as.character(1:nrow(dlbcl.data)),sep="")

# create train and test data objects. censoring.status=1 means the event occurred;
#  censoring.status=0 means censored

data<-list(x=dlbcl.train.data,y=dlbcl.train.time, censoring.status=dlbcl.train.cens, featurenames=featurenames)
data.test<-list(x=dlbcl.test.data,y=dlbcl.test.time, censoring.status=dlbcl.test.cens, featurenames= featurenames)


# train  the model. This step just computes the scores for each feature
train.obj<- superpc.train(data, type="survival")

# note for regression (non-survival) data, we leave the component "censoring.status"
# out of the data object, and call superpc.train with type="regression".
# otherwise the superpc commands are all the same


# cross-validate the model (i.e., identify Cox score threshold at which markers should be retained)
cv.obj<-superpc.cv(train.obj, data)

#plot the cross-validation curves. From this plot we see that the 1st 
# principal component is significant and the best threshold  is around 1.5
superpc.plotcv(cv.obj)


# here we have the luxury of  test data, so we can compute the likelihood ratio statistic
# over the test data and plot them. We see that the threshold of 1.5
# works pretty well (see JASA paper, e.g., section 5)
lrtest.obj<-superpc.lrtest.curv(train.obj, data, data.test)
superpc.plot.lrtest(lrtest.obj)


# now we derive the predictor of survival for the test data, 
# and then then use it
# as the predictor in a Cox model . We see that the 1st supervised PC is
# highly significant; the last 2 are not
fit.cts<- superpc.predict(train.obj, data, data.test, threshold=1.5, n.components=3, prediction.type="continuous")
superpc.fit.to.outcome(train.obj, data.test, fit.cts$v.pred)


# sometimes a discrete (categorical) predictor is attractive.
# Here we form two groups by cutting the predictor at its median
#and then plot Kaplan-Meier curves for the two groups
fit.groups<- superpc.predict(train.obj, data, data.test, threshold=1.0, n.components=1, prediction.type="discrete")
superpc.fit.to.outcome(train.obj, data.test, fit.groups$v.pred)

plot(survfit(Surv(data.test$y,data.test$censoring.status)~fit.groups$v.pred), col=2:3, xlab="time", ylab="Prob survival")


# Finally, we look for a predictor of survival a small number of
# markers (rather than all p markers). We do this by computing an importance
# score (see JASA paper, section 4) for each marker
# (if data standardized, equal to its correlation with the first supervised PC predictor).
# Features with the highest correlation contribute
# most to the prediction of the outcome.
# We soft threshold the importance scores, and use the shrunken
# scores as marker weights to form a reduced predictor. Cross-validation
# gives us an estimate of the best amount to shrink and an idea of
# how well the shrunken predictor works.
fit.red<- superpc.predict.red(train.obj, data, data.test, threshold=1.5)
fit.redcv<- superpc.predict.red.cv(fit.red, cv.obj,  data,  threshold=1.5)
superpc.plotred.lrtest(fit.redcv)


# Finally we list the significant markers, in order of decreasing importance score
# raw-score is marginal Cox score
superpc.listfeatures(data.test, train.obj, fit.red)



#A note on interpretation:
#The signs of the scores (latent factors) v.pred returned by superpc.predict are chosen so that the regression of the outcome on each factor has a positive coefficient. This is just a convention to aid in interpretation.
#
#
#For survival data, this means
#Higher score => higher risk (worse survival)
#
#For regression data, this means
#Higher score => higher mean of the outcome
#
#How about the direction of effect for each individual feature (gene)? The function superpc.listfeatures reports an importance score equal to the correlation between each feature and the latent factor.
#
#Hence for survival data,
#
#Importance score positive means
#increase in value of feature=> higher risk (worse survival)
#
#For regression data,
#
#Importance score positive means
#increase in value of feature => higher mean of the outcome
#
#The reverse are true for Importance score negative. 
#
