#Random forest - both data types (phylo eigenvetors and ecomorph traits), few variables
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

#select cols
kp <-c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:5, sep=""), paste("source", "eig", 1:5, sep="")) #cut to 21 because there are 21 ecomorphological variables
nms <- c("sourceTaxonName","targetTaxonName","targetBodyMass.Value","sourceBodyMass.Value","sourceeat_plants","sourceeat_invs","targetflight","sourceForStrat.wataroundsurf","targeteat_plants","targeteat_verts","targetActivity.Crepuscular","target_aerial_mam")
kp <- unique(c(kp,nms))
Int <- Int[,names(Int)%in%c(kp,"interact","outside")]
Non <- Non[,names(Non)%in%c(kp,"interact","outside")]
SD_foc <- SD_foc[,names(SD_foc)%in%c(kp,"interact","outside")]

df_opt_thrw <- replicate(100,rep_all_thrw(Int, Non, Non,ins=0.75,out=1,abs=2.5,thresh=0.41,mtr=3,max.depth=1000, trees=800)) #mean TSS~0.758 (if use 10 eig, 0.731; 20 = 0.743; 6=0.723)
saveRDS(df_opt_thrw, "###/GloBIonly_bothDataTypes_fewVar.rds")

#########################################################################################################################
SDall <- replicate(100,rf_thresh_all(Int,Non,ins=0.75,out=1,abs=2.5,SD_foc,thresh=0.41,mtr=3,max.depth=1000,num.trees=800)) #
mean(unlist(SDall[1,])) #TSS 
mean(unlist(SDall[2,])) #AUC 
saveRDS(SDall, "###/SDaspplied_bothDataTypes_fewVar.rds")

