#randomization tests comparing random forest performance predicting predator-prey interactions for native versus recently introduced (cats and foxes) predators)
#set file paths to run (search for lines with "###")

library(simstudy)
library(data.table)

setwd("~/###")
nats <-readRDS("###/SDaspplied_ecomorph_nativePreds.rds")
cf <- readRDS("###/SDaspplied_ecomorph_manyVar_cats&foxes.rds")

#reformat SD data
nats <- unlist(nats[1,])
cf <- unlist(cf[1,])

#from https://alfurka.github.io/2020-12-16-sampling-with-density-function/
samp_func <- function (x, N.draw){
dens = density(x)
cdf = cumsum(dens$y) / cumsum(dens$y)[length(dens$y)] #a random draw from the cumulative distribution function (CDF) that finds the index of the draw. I add 1 to the index because it finds the index value on the left. From ##dens$x[] I find the value of the drawn index.
samp = replicate(N.draw, dens$x[findInterval(runif(1), cdf)+1])
return(samp)}

#Comparing models with different variables
#GloBI trained and tested
natsSamp <- samp_func(nats,100000)
cfSamp <- samp_func(cf,100000)

dt <- data.frame(cbind(natsSamp,cfSamp))
dt$natsVScf <- dt$natsSamp>=dt$cfSamp
dt$cfVSnats <- dt$cfSamp>=dt$natsSamp

table(dt$natsVScf)
table(dt$cfVSnats)

#get probability:
table(dt$cfVSnats)[1]/100000
