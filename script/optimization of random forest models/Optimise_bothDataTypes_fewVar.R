#Optimise random forest - model with few variables (both phylogenetic eigenvectors and ecomorph traits)
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

#select columns to keep
kp <-c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:5, sep=""), paste("source", "eig", 1:5, sep="")) #cut to 21 because there are 21 ecomorphological variables
nms <- names(Non)[!(grepl("eig",names(Non)))]
kp <- unique(c(kp,nms))

Int <- Int[,names(Int)%in%kp]
Non <- Non[,names(Non)%in%kp]

#check optimal number of unobserved (abs) versus observed interactions 
variables = 2
iterations = 20
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:20) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=0.25,out=1,abs=i/2,thresh=0.4,mtr=14,max.depth=1000,num.trees = 800)))
  output[i,] <- c(round(dt,3),i/2)
  print(dt)
}

output=data.frame(output)
pdf("###/OpPEMs+Eco_abs_fewVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) #peaks at 
dev.off()

absv <- output$X2[which.max(output$X1)]

#check how many absent observations from inside versus outside suitable range is best. 
variables = 2
iterations = 10
output <- matrix(ncol=variables,nrow=iterations)
for(i in 1:10) {
  dt <- mean(replicate(40,rep_all_thrw(Int,Non,Non,ins=i/4,out=1,abs=absv,thresh=0.42,mtr=6,max.depth=1000,num.trees=800)))
  output[i,] <- c(round(dt,3),i/4)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs+Eco_ins-outs_fewVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) #peaks at
dev.off()

insv <- output$X2[which.max(output$X1)]

#optimise mtry
variables = 2
iterations = 20
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:20) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int, Non, Non,ins=insv,out=1,abs=absv,thresh=0.42,mtr=i,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),i)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs+Eco_mtry_fewVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) 
dev.off()

mtryv <- output$X2[which.max(output$X1)]

#optimise thresh
variables = 2
iterations = 15
output <- matrix(ncol=variables,nrow=iterations)
for (i in 1:15) {
  dt <- mean(replicate(40,rep_all_thrw_trees(Int,Non,Non,ins=insv,out=1,abs=absv,thresh=0.36+i/100,mtr=mtryv,max.depth=1000,trees = 800)))
  output[i,] <- c(round(dt,3),0.36+i/100)
  print(i)
}
output=data.frame(output)
pdf("###/OpPEMs_thresh_fewVar.pdf", height = 8, width = 16)
plot(output$X2,output$X1) #peaks at 0.42
dev.off()

thrv <- output$X2[which.max(output$X1)]
