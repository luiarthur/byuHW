gradespre <- read.table("grade.dat", col.names=
   c("Labmean","Labquiz","Hwmean","Pqmean","Project","Exam1","Exam2","ExamFin"))
grades <- gradespre[,-c(2,5)]
n <- nrow(grades)
p <- ncol(grades)
grades

X <- grades

## PCA on original data ##
eiginfo.X <- eigen(var(X))
eiginfo.X
## % and cumulative % for each principal component ##
cbind(eiginfo.X$values/sum(eiginfo.X$values),
      cumsum(eiginfo.X$values)/sum(eiginfo.X$values))




##### Selecting number of PCs (using original data) #####

### Approach (1) ###
cbind(eiginfo.X$values/sum(eiginfo.X$values),
      cumsum(eiginfo.X$values)/sum(eiginfo.X$values))
### Approach (1) retains 2 components for k=.80 and 2 components for k=.90

### Approach (2) ###
mean(eiginfo.X$values)
eiginfo.X$values
### Approach (2) retains 1 components

### Approach (3) ###
ts.plot(eiginfo.X$values,xlab="component",ylab="e'val")
abline(h=mean(eiginfo.X$values))
title("Scree Plot")
### Approach (3) retains 1? component 



##############################
##############################
## PCA on standardized data ##
##############################
##############################
eiginfo.X <- eigen(cor(X))
eiginfo.X
## % and cumulative % for each principal component ##
cbind(eiginfo.X$values/sum(eiginfo.X$values),
      cumsum(eiginfo.X$values)/sum(eiginfo.X$values))



##### Selecting number of PCs (using original data) #####

### Approach (1) ###
cbind(eiginfo.X$values/sum(eiginfo.X$values),
      cumsum(eiginfo.X$values)/sum(eiginfo.X$values))
### Approach (1) retains 3 components for k=.80 and 5 components for k=.90

### Approach (2) ###
mean(eiginfo.X$values)
eiginfo.X$values
### Approach (2) retains 2 components

### Approach (3) ###
ts.plot(eiginfo.X$values,xlab="component",ylab="e'val")
abline(h=mean(eiginfo.X$values))
title("Scree Plot")
### Approach (3) retains 1? component 
