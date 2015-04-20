#Attached is the data set of the scraped reviews.  It includes two data frames:
#train and test.  Use the data in train to create your model.  Use the data in
#test to get a test mean squared error.  
#
#Rules: 
#1. Use the content of the data however you want (fit splines, count #words, etc.)
#2. You must fit a model using the odds/probabilities -- nominal or
#   ordinal -- we discussed in class (but you can use different links, just be sure
#   you know what the model is you're fitting) 
#4. Note that the data size is actually quite small for fitting four
#   different models, so be lenient with level of significance.  
#5. Write a 1-2 page report including your model, interpretation of parameters, 
#   prediction error, plots of interest, and a brief review of why this model 
#   is valuable.  
#6. In class, each person will present their model so that we can vote on the
#   winners for the various categories. I'd like this to be brief, so no slides
#   are needed unless necessary (I'll let you decide if that's the case).
#7. There will be four categories for prizes: 
#     - Smallest test error predictions (in terms of mean squared error)
#     - Model that tells the best story (quantifies relationships in a valuable/interesting way)
#     - Most interpretable (but useful) model
#     - Most creative model (e.g., uses statistical words in the model, etc.).
#   The class will vote on the last three categories.  



# Functions:
  rgx <- function(pattern,text) {
    obj <- gregexpr(pattern,text)[[1]] 
    pos <- as.numeric(obj)
    counts <- length(pos)
    
    out <- list("counts"=counts,"position"=pos)
    out
  }

  wd.count <- function(text) sapply(gregexpr("\\W+", text), length) + 1

# Read Data:
  load("reviews.Rdata") # test, train
  train.review <- sapply(train$review,as.character)

  pairs(train[,-8])
  X <- train[,2:6]
  y <- train[,1]
  pairs(X)

# Explore:
  star.summary <- function(k,stars=train$star,reviews=train.review) {
    ind <- which(stars==k)
    words <- paste(reviews[ind],collapse=" ")
    word.tab <- table(strsplit(words," "))
    word.tab <- word.tab[order(word.tab,decreasing=TRUE)]

    list("tab"=word.tab,"counts"=wd.count(reviews[ind]))
  }

  star.sum.1 <- star.summary(1)
  star.sum.2 <- star.summary(2)
  star.sum.3 <- star.summary(3)
  star.sum.4 <- star.summary(4)
  star.sum.5 <- star.summary(5)

  star.sum.1$counts
  star.sum.2$counts
  star.sum.3$counts
  star.sum.4$counts
  star.sum.5$counts

  star.sum.1$tab[star.sum.1$tab>=3]
  star.sum.2$tab[star.sum.2$tab>=3]
  star.sum.3$tab[star.sum.3$tab>=3]
  star.sum.4$tab[star.sum.4$tab>=3]
  star.sum.5$tab[star.sum.5$tab>=3]

  neg.words <- "too|nothing|why|how|shortcoming|didn't like|don't like|difficult|stinks|sucks|blows|annoy|bad|not recommend"
  pos.words <- "excellent|everything|easy to follow|very good|mathematical|intuition|d recommend this|theoretical"

  neg.words.count <- apply(as.matrix(train.review),1,function(x) rgx(neg.words,x)$counts)
  pos.words.count <- apply(as.matrix(train.review),1,function(x) rgx(pos.words,x)$counts)

  X <- cbind(X,neg.words.count,pos.words.count)
  X <- cbind(X[,1]/X[,2],X)
  colnames(X) <- c("helpful.prop","helpful","help.count","day","month","year","neg.count","pos.count")
# Fit model:
  library(VGAM)
  pairs(cbind(y,X))
  x <- as.data.frame(X)
  head(X)

  #mod <- vglm(y~helpful.prop+helpful+help.count+month+year+bs(neg.count)+bs(pos.count),family=cumulative(parallel=T),data=x)
  #mod <- vglm(y~helpful.prop+helpful+help.count+month+year+bs(neg.count),family=cumulative(parallel=T),data=x)
  #mod <- vglm(y~helpful.prop+helpful+help.count+year+bs(neg.count),family=cumulative(parallel=T),data=x)
  mod <- vglm(y~bs(helpful.prop)+year+bs(neg.count),family=cumulative(parallel=T),data=x)
  summary(mod)
  
  pairs(cbind(jitter(y),x$helpful.prop,x$year,jitter(x$neg.count)),
        labels=c("y","help prop","year","neg count"),
        main="Paired Scatter Plots",col="blue")

# Test Model:
  create.X <- function(m) {
    #m <- test[,-1]
    neg.wd.count <- apply(as.matrix(m$review),1,function(x) rgx(neg.words,x)$counts)
    pos.wd.count <- apply(as.matrix(m$review),1,function(x) rgx(pos.words,x)$counts)

    m <- data.frame(m$help.yes/m$help.total,m$year,neg.wd.count)
    colnames(m) <- c("helpful.prop","year","neg.count")

    #m <- data.frame(m[,1]/m[,2],m[,1],m[,2],as.factor(m[,4]),m[,5],neg.wd.count,pos.wd.count)
    #colnames(m) <- c("helpful.prop","helpful","help.count","month","year","neg.count","pos.count")

    #m <- data.frame(m[,1]/m[,2],m[,1],m[,2],as.factor(m[,4]),m[,5],neg.wd.count)
    #colnames(m) <- c("helpful.prop","helpful","help.count","month","year","neg.count")
    
    #m <- data.frame(m[,1]/m[,2],m[,1],m[,2],m[,5],neg.wd.count)
    #colnames(m) <- c("helpful.prop","helpful","help.count","year","neg.count")

    m
  }
  
  

  yi <- test[-16,1] # 16 was missing
  pred <- predict(mod,create.X(test[,-1]),type="response")
  y.hat <- apply(pred,1,function(x) which.max(x))
  
  table(yi,y.hat)
  cbind(yi,y.hat)

  (mse <- mean((yi-y.hat)^2))
