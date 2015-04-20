source("genDat2.R")
source("gibbs.R")
source("ibp.R")

elapsed.time <- system.time(out <- gibbs.post(Y,a=1,B=2000,burn=0,showProgress=T,
                                              plotProgress=T,a.a=3,a.b=2,
                                              siga=1,sigx=.5))

