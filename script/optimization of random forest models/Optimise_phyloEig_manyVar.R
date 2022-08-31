#Optimise random forest - model with many variables (phylogenetic eigenvectors only)
#To start with, set parameters to ins=1, thresh=0.5, mtr=round(sqrt(Vars)). Update values between loops below, and run through loops twice to optimise parameters. 
#set file paths to run (search for lines with "###")

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)
library(caTools)

Int <- readRDS("###/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/allNon_sameCont.RDS")

#select columns to keep
kp <-c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:21, sep=""), paste("source", "eig", 1:21, sep="")) #cut to 21 because there are 21 ecomorphological variables
Int <- Int[,names(Int)%in%kp]
Non <- Non[,names(Non)%in%kp]


#check optimal number of unobserved (abs) versus observed interactions 
variables = 2
iterations = 20
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:20) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=0.75,out=1,abs=i/2,thresh=0.33,mtr=1,max.depth=1000,num.trees = 800)))
  output[i,] <- c(round(dt,3),i/2)
  print(i)
}

output=data.frame(output)
pdf("###/OpPEMs_abs_manyVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

absv <- output$X2[which.max(output$X1)]

#check how many absent observations from inside versus outside suitable range is best. 
variables = 2
iterations = 10
output <- matrix(ncol=variables,nrow=iterations)
for(i in 1:10) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=i/4,out=1,abs=absv,thresh=0.35,mtr=1,max.depth=1000,num.trees=800)))
  output[i,] <- c(round(dt,3),i/4)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs_ins-outs_manyVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) #peaks at 1:2
dev.off()

insv <- output$X2[which.max(output$X1)]


#optimise mtry
variables = 2
iterations = 10
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:10) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int, Non, Non,ins=insv,out=1,abs=absv,thresh=0.35,mtr=i,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),i)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs_mtry_manyVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

mtryv <- output$X2[which.max(output$X1)]

#optimise thresh
variables = 2
iterations = 15
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:15) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int,Non,Non,ins=insv,out=1,abs=absv,thresh=0.29+i/100,mtr=mtryv,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),0.29+i/100)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs_thresh_manyVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

thrv <- output$X2[which.max(output$X1)]

