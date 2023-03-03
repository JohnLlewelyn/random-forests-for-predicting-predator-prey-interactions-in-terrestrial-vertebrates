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

gdata <- RF_allInt(x,y,ins,out,abs)
