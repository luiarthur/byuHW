source("rfunctions.R")
source("ibp.R")
source("generateData.R")
source("gibbs.R")

elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=1000,burn=0,showProgress=T,
                                              plotProgress=T,a.a=1,a.b=1))

# What Next?

