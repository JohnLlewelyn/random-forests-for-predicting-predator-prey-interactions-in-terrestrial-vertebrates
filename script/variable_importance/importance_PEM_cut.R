#identify important variables
# "###" indicates where filepaths are needed
source("###/opt_functions.R")
source("###/variable_importance_functions.R")

#libraries
library(ranger)
library(caTools)
library(tidyverse)

#get  global dataset 
Int <- readRDS("###/data/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/data/allNon_sameCont.RDS")

#cut to species with 5 or more records
ch<-data.frame(table(Int$sourceTaxonName))
prds <- ch$Var1[ch$Freq>4]
Int<-Int[Int$sourceTaxonName%in%prds,]
Non<-Non[Non$sourceTaxonName%in%prds,]

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0

#cut to required columns
kp <- c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:21, sep=""), paste("source", "eig", 1:21, sep="")) #cut to 21 because there are 21 ecomorphological variables
Int <- Int[,names(Int)%in%c(kp,"interact","outside")]
Non <- Non[,names(Non)%in%c(kp,"interact","outside")]

#set variable for building dataset
x=Int
y=Non
abs=4.75
ins=1
out=1
#build dataset
set.seed(456)
gdata <- RF_allInt(x,y,ins,out,abs)

#fit RF and get importance
dat=gdata[[1]]
wgt=gdata[[2]]

#run a RF and get the least important variable, remove from data and rerun until 10 variables remain (remove 32 variables in this case)
#1st
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#2st
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#3rd
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#4th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#5th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#6th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#7th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#8th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#9th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#10th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#11th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#12th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#13th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#14th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#15th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]
#16th
set.seed(456)
vrem <-least_imp_default(dat, wgt)
vrem
vrem <- c(paste("target",vrem, sep=""), paste("source",vrem, sep=""))
dat <- dat[,-which(names(dat)%in%vrem)]

names(dat)
saveRDS(dat,"###/keepVar_PEM_cut.RDS")

