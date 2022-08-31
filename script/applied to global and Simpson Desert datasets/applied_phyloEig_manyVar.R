#Random forest - using many (21) phylogenetic eignevectors. Apply 100 times to GloBI and Simpson Desert data
#set file paths to run (search for lines with "###")

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)
library(caTools)

Int <- readRDS("###/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/allNon_sameCont.RDS")
SD_foc <- readRDS("###/allperms_cut2_20EVs.RDS")

#select required columns
nms <-c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:21, sep=""), paste("source", "eig", 1:21, sep="")) #
Int <- Int[,names(Int)%in%c(nms,"interact","outside")]
Non <- Non[,names(Non)%in%c(nms,"interact","outside")]

df_opt_thrw <- replicate(100,rep_all_thrw(Int, Non, Non,ins=0.75,out=1,abs=6.5,thresh=0.33,mtr=1,max.depth=1000, trees=800)) #mean TSS~0.758 (if use 10 eig, 0.731; 20 = 0.743; 6=0.723)
saveRDS(df_opt_thrw, "###/GloBIonly_phyloEig_manyVar.rds")

#########################################################################################################################
SDall <- replicate(100,rf_thresh_all(Int,Non,ins=0.75,out=1,abs=6.5,SD_foc,thresh=0.33,mtr=1,max.depth=1000,num.trees=800)) #
mean(unlist(SDall[1,])) #TSS 
mean(unlist(SDall[2,])) #AUC
saveRDS(SDall, "###/SDaspplied_phyloEig_manyVar.rds")

