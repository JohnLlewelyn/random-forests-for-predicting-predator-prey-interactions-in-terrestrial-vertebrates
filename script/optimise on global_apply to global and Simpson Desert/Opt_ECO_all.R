#optimise using ecomorph,many variables, and the full/original global dataset 
# "###" indicates where filepaths need to be added
source("###/opt_functions.R")

#libraries
library(ranger)
library(caTools)
library(tidyverse)

Int <- readRDS("###/data/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/data/allNon_sameCont.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0

#remove eig 
nms <- names(Non)[!(grepl("eig",names(Non)))]
Int <- Int[,names(Int)%in%c(nms,"interact","outside")]
Non <- Non[,names(Non)%in%c(nms,"interact","outside")]

thresh=0.5
n_features=length(names(Int[,-tax_cols(Int)]))
mtr=round(sqrt(n_features))
ntrees = 100
mdepth = 0
ins=1
out=1


#check optimal number of unobserved (abs) versus observed interactions 
set.seed(123)
a <- lapply(seq(0.25,8,0.25), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=x,ins=ins,out=out,thresh=thresh,mtr=mtr,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(seq(0.25,8,0.25), unlist(a))
plot(output[,1],output[,2],main="abs")
absv = output[,1][which.max(output[,2])]
output_abs1=output
print(absv)


#check how many absent observations from inside versus outside suitable range is best. 
set.seed(123)
a <- lapply(seq(0.25,4,0.25), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=x,out=out,thresh=thresh,mtr=mtr,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(seq(0.25,4,0.25), unlist(a))
plot(output[,1],output[,2],main="ins:out")
insv = output[,1][which.max(output[,2])]
output_ins1=output
print(insv)

#optimise mtry
set.seed(123)
a <- lapply(2:round(n_features*.6), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=thresh,mtr=x,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(2:round(n_features*.6), unlist(a))
plot(output[,1],output[,2],main="mtry")
mtryv = output[,1][which.max(output[,2])]
output_mtry1=output
print(mtryv)

#optimise thresh
set.seed(123)
a <- lapply(seq(0.3,0.6,0.01), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=x,mtr=mtryv,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(seq(0.3,0.6,0.01), unlist(a))
plot(output[,1],output[,2],main="thresh")
threshv = output[,1][which.max(output[,2])]
output_thresh1=output
print(threshv)

#optimise num.trees
set.seed(123)
a <- lapply(seq(50,800,50), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=threshv,mtr=mtryv,mdepth = mdepth, ntrees = x))))})

output <- cbind(seq(50,800,50), unlist(a))
plot(output[,1],output[,2],main="ntrees")
ntreesv = output[,1][which.max(output[,2])]
output_trees1=output
print(ntreesv)

#optimise tree depth
set.seed(123)
a <- lapply(c(10,20,30,50,100,500,1000,0), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=threshv,mtr=mtryv,mdepth = x, ntrees = ntreesv))))})

output <- cbind(c(10,20,30,50,100,500,1000,0), unlist(a))
plot(output[,1],output[,2],main="max depth")
mdepthv = output[,1][which.max(output[,2])] #
output_depth1=output
print(mdepthv)

####re-try to tweak############################################################################################
#updated parameters
thresh=threshv
mtr=mtryv
ntrees = ntreesv
mdepth = mdepthv
ins=insv
out=1
#abs=absv# don't need to set this because it's the first optimised paramter
#check optimal number of unobserved (abs) versus observed interactions 

rm(list = c('mtryv','ntreesv','mdepthv','insv','threshv'))

#abs
set.seed(123)
a <- lapply(seq(0.25,8,0.25), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=x,ins=ins,out=out,thresh=thresh,mtr=mtr,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(seq(0.25,8,0.25), unlist(a))
plot(output[,1],output[,2],main="abs")
absv = output[,1][which.max(output[,2])]
output_abs2=output
print(absv)

#check how many absent observations from inside versus outside suitable range is best. 
set.seed(123)
a <- lapply(seq(0.25,4,0.25), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=x,out=out,thresh=thresh,mtr=mtr,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(seq(0.25,4,0.25), unlist(a))
plot(output[,1],output[,2],main="ins:out")
insv = output[,1][which.max(output[,2])]
output_ins2=output
print(insv)

#optimise mtry
set.seed(123)
a <- lapply(2:round(n_features*.6), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=thresh,mtr=x,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(2:round(n_features*.6), unlist(a))
plot(output[,1],output[,2],main="mtry")
mtryv = output[,1][which.max(output[,2])]
output_mtry2=output
print(mtryv)

#optimise thresh
set.seed(123)
a <- lapply(seq(0.3,0.6,0.01), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=x,mtr=mtryv,mdepth = mdepth, ntrees = ntrees))))})

output <- cbind(seq(0.3,0.6,0.01), unlist(a))
plot(output[,1],output[,2],main="thresh")
threshv = output[,1][which.max(output[,2])]
output_thresh2=output
print(threshv)

#optimise num.trees
set.seed(123)
a <- lapply(seq(50,800,50), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=threshv,mtr=mtryv,mdepth = mdepth, ntrees = x))))})

output <- cbind(seq(50,800,50), unlist(a))
plot(output[,1],output[,2],main="ntrees")
ntreesv = output[,1][which.max(output[,2])]
output_trees2=output
print(ntreesv)

#optimise tree depth
set.seed(123)
a <- lapply(c(10,20,30,50,100,500,1000,0), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=threshv,mtr=mtryv,mdepth = x, ntrees = ntreesv))))})

output <- cbind(c(10,20,30,50,100,500,1000,0), unlist(a))
plot(output[,1],output[,2],main="max depth")
mdepthv = output[,1][which.max(output[,2])] #
output_depth2=output
print(mdepthv)

set.seed(123)
eco_perf_global <- replicate(100,optRF_allM_depth(Int,Non,Non,abs=absv,ins=insv,out=out,thresh=threshv,mtr=mtryv,mdepth = mdepthv, ntrees = ntreesv))
glb = data.frame(t(reduce(eco_perf_global, full_join, by = "X1")))
names(glb) <- glb[1,]
glb<-glb[-1,]
glb <- sapply(glb, as.numeric) 
glb <- as.data.frame(glb) #
#mean tss


#save results
saveRDS(glb,"###/results/EcoAll_GLB_rp5_100.RDS")

#now train on all the global data and apply to Simpson Desert#######
SD_foc <- readRDS("###/data/allperms_cut2_20EVs.RDS")
SD_foc$source_aerial_mam <- 0

set.seed(123)
SD_perf <- replicate(100, RF_SD_simp_depth(Int,Non,dt_test=SD_foc, ins=insv, out=out,abs=absv,thresh=threshv,mtr=mtryv,mdepth=mdepthv,ntrees=ntreesv))
SDp = data.frame(t(reduce(SD_perf, full_join, by = "X1")))
names(SDp) <- SDp[1,]
SDp<-as.data.frame(SDp[-1,])
SDp <- as.data.frame(sapply(SDp, as.numeric)) 
#mean tss

saveRDS(SDp,"###/results/EcoAll_SD_rp5_100.RDS")

#save hyper parameters
hps <- data.frame(Hpar=c("abs","ins","mtry","thresh","ntrees","depth"),val=c(absv,insv,mtryv,threshv,ntreesv,mdepthv))
saveRDS(hps,"###/results/Hpars_EcoAll.RDS")


#script for figure to show why ins:outs is important
set.seed(123)
a <- lapply(seq(0.25,4,0.25), function(x){
  (mean(replicate(5,optRF_depth(Int,Non,Non,abs=absv,ins=x,out=out,thresh=threshv,mtr=mtryv,mdepth = mdepthv, ntrees = ntreesv))))})

output <- cbind(seq(0.25,4,0.25), unlist(a))

plot(output[,1],output[,2],pch=19,ylab="TSS",xlab="ratio of non-interactions from\ninside:outside preferred size range")
