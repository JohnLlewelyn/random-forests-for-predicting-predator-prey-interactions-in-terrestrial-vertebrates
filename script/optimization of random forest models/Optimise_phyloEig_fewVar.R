#Optimise random forest - model with few variables (phylogenetic eigenvectors only)
#To start with, set parameters to ins=1, thresh=0.5, mtr=round(sqrt(Vars)). Update values between loops below, and run through loops twice to optimise parameters. 
#set file paths to run (search for lines with "###")

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)
library(caTools)

Int <- readRDS("###/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/allNon_sameCont.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0

#select columns to keep
nms <-c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:5, sep=""), paste("source", "eig", 1:5, sep="")) #
Int <- Int[,names(Int)%in%c(nms,"interact","outside")]
Non <- Non[,names(Non)%in%c(nms,"interact","outside")]

#check optimal number of unobserved (abs) versus observed interactions 
variables = 2
iterations = 15
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:15) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=3,out=4,abs=i/2,thresh=0.46,mtr=7,max.depth=1000,num.trees = 800)))
  output[i,] <- c(round(dt,3),i/2)
  print(i)
}

output=data.frame(output)
pdf("###/OpPEMs21_abs_fewVars.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

absv <- output$X2[which.max(output$X1)]

#check how many absent observations from inside versus outside suitable range is best. 
variables = 2
iterations = 10
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:10) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=i/4,out=1,abs=absv,thresh=0.46,mtr=7,max.depth=1000,num.trees=800)))
  output[i,] <- c(round(dt,3),i/4)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs21_ins-outs_fewVars.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

insv <- output$X2[which.max(output$X1)]

#optimise mtry
variables = 2
iterations = 9
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:9) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int, Non, Non,ins=insv,out=1,abs=absv,thresh=0.46,mtr=i,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),i)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs21_mtry_fewVars.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

mtryv <- output$X2[which.max(output$X1)]


#optimise thresh
variables = 2
iterations = 15
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:15) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int,Non,Non,ins=insv,out=1,abs=absv,thresh=0.35+i/100,mtr=mtryv,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),0.35+i/100)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs5_thresh_fewVars.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

