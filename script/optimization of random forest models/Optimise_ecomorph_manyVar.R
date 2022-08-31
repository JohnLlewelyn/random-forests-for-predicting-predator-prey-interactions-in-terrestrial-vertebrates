#Optimise random forest - model with many variables (ecomorph traits only)
#To start with, set parameters to ins=1, thresh=0.5, mtr=round(sqrt(Vars)). Update values between loops below, and run through loops twice to optimise parameters. 
#set file paths to run (search for lines with "###")

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)
library(caTools)

Int <- readRDS("data/GloBIplus_Int20EVs.RDS")
Non <- readRDS("data/allNon_sameCont.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0

#remove eig > 600
nms <- names(Non)[!(grepl("eig",names(Non)))]
Int <- Int[,names(Int)%in%c(nms,"interact","outside")]
Non <- Non[,names(Non)%in%c(nms,"interact","outside")]

#check optimal number of unobserved (abs) versus observed interactions 
variables = 2
iterations = 25
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:25) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=0.5,out=1,abs=i/2,thresh=0.42,mtr=6,max.depth=1000,num.trees = 800)))
  output[i,] <- c(round(dt,3),i/2)
  print(dt)
}

output=data.frame(output)
absv <- output$X2[which.max(output$X1)]

pdf("###/OpECO_abs.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

#check how many absent observations from inside versus outside suitable range is best. 
variables = 2
iterations = 10
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:10) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=i/4,out=1,abs=absv,thresh=0.42,mtr=12,max.depth=1000,num.trees=800)))
  output[i,] <- c(round(dt,3),i/4)
  print(dt)
}
output=data.frame(output)

insv <- output$X2[which.max(output$X1)]

pdf("###/OpECO_ins-outs.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

#optimise mtry
variables = 2
iterations = 20
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:20) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int, Non, Non,ins=insv,out=1,abs=absv,thresh=0.42,mtr=i,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),i)
  print(dt)
}
output=data.frame(output)

mtryv <- output$X2[which.max(output$X1)]

pdf("###/OpECO_mtry.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

#optimise thresh
variables = 2
iterations = 11
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:11)  {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int,Non,Non,ins=insv,out=1,abs=absv,thresh=0.40+i/100,mtr=mtryv,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),0.40+i/100)
  print(dt)
}
output=data.frame(output)
pdf("###/OpECO_thresh.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()
