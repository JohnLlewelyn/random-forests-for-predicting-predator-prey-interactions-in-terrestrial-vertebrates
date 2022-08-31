#Random forest - many ecomorph traits. Apply 100 times to GloBI and Simpson Desert data
#set file paths to run (search for lines with "###") 

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)
library(caTools)

Int <- readRDS("###/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/allNon_sameCont.RDS")
SD_foc <- readRDS("###/allperms_cut2_20EVs.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0
SD_foc$source_aerial_mam <- 0

#remove cols
kp <- names(Non)[!(grepl("eig",names(Non)))]
Int <- Int[,names(Int)%in%c(kp,"interact","outside")]
Non <- Non[,names(Non)%in%c(kp,"interact","outside")]
SD_foc <- SD_foc[,names(SD_foc)%in%c(kp,"interact","outside")]

df_opt_thrw <- replicate(100,rep_all_thrw(Int, Non, Non,ins=1,out=1,abs=4,thresh=0.43,mtr=11,max.depth=1000, trees=800)) 
saveRDS(df_opt_thrw, "###/GloBIonly_ecomorph_manyVar.rds")

#########################################################################################################################
SDall <- replicate(100,rf_thresh_all(Int,Non,ins=1,out=1,abs=4,SD_foc,thresh=0.43,mtr=11,max.depth=1000,num.trees=800)) #
mean(unlist(SDall[1,])) #TSS 
mean(unlist(SDall[2,])) #AUC 
saveRDS(SDall, "###/SDaspplied_ecomorph_manyVar.rds")

