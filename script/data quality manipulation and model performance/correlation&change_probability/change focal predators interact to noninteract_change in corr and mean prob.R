#Training dataset taxonomic coverage (TDTC) and performance - change predator interactions to non-interactions 
#(simulating cases where predator interactions understudied and high prevalence of false negatives)
#set file paths (on lines with "###")

setwd("###")
source("###/all_functions_ranger.R")
source("###/opt_functions.R")

#libraries
library(ranger)
library(plyr)

#get data
Int <- readRDS("###/data/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/data/allNon_sameCont.RDS")
SD_foc <- readRDS("###/data/allperms_cut2_20EVs.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0
SD_foc$source_aerial_mam <- 0

#cut global dataset to species with 5 or more records
ch<-data.frame(table(Int$sourceTaxonName))
prds <- ch$Var1[ch$Freq>4]
Int<-Int[Int$sourceTaxonName%in%prds,]
Non<-Non[Non$sourceTaxonName%in%prds,]

#get lists of Simpson Desert predators, prey (which has both)
prds <- unique(SD_foc$sourceTaxonName)
prey <- unique(SD_foc$targetTaxonName)

#remove columns not required
kp <- c("targetTaxonName","sourceTaxonName","interact","outside",paste("target", "eig", 1:21, sep=""), paste("source", "eig", 1:21, sep="")) #cut to 21 because there are 21 ecomorphological variables
nms <- names(Non)[!(grepl("eig",names(Non)))]
kp <- unique(c(kp,nms))
Int <- Int[,names(Int)%in%c(kp,"interact","outside")]
Non <- Non[,names(Non)%in%c(kp,"interact","outside")]
SD_foc <- SD_foc[,names(SD_foc)%in%c(kp)]  #deleted "interact","outside" because don't need this in this dataset

#function so can get correlations split by a grouping column (predator)
split_corr <- function(xx)
{
  return(data.frame(COR = cor(xx$predict_orig, xx$predict_reduced)))
}

#function so can get change in mean probability for each predator
split_mean <- function(xx)
{
  return(data.frame(mean_diff=-mean(xx$predict_orig)+mean(xx$predict_reduced)))
}

#the function for changing predators' interact to noninteract
prdChange_TDTC_corrMOD <- function(x,y,ins,out,abs,dt_test,mtr,ntrees,mdepth,percPred,...){
  obs <-  x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$obs = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  data$outside = NULL
  data_orig <- data
  data_orig <- data_orig[,-tax_cols(data_orig)]
  data_w1 <- ifelse(data_orig$interact=="TRUE",1,1/(table(data_orig$interact)[2]/table(data_orig$interact)[1]))
  #split training data so SD_preds in one group
  SDprd <- data[data$sourceTaxonName%in%prds&data$interact==TRUE,]
  Oprd <- data[!(data$sourceTaxonName%in%prds&data$interact==TRUE),]
  #remove rows based on perc (remove from data and data_w)
  sel1 <- round(percPred/100*length(SDprd$sourceTaxonName))
  sel2 <- sample(1:length(SDprd$sourceTaxonName),sel1)
  SDprd$interact[] <- 'FALSE'
  SDprd$interact[sel2] <- 'TRUE'
  data2 <- rbind(SDprd,Oprd)
  data_w2 <- ifelse(data2$interact=="TRUE",1,1/(table(data2$interact)[2]/table(data2$interact)[1]))
  dt_test = dt_test[,names(data2)]
  dt_test$interact <- as.factor(dt_test$interact)
  data2 = data2[,-tax_cols(data2)]
  rf1 = ranger(formula = interact ~., 
               data = data_orig,
               mtry = mtr, 
               num.threads = 20, 
               probability = T, 
               importance = 'impurity',
               case.weights = data_w1,
               num.trees = ntrees,
               max.depth = mdepth)
  rf2 = ranger(formula = interact ~., 
               data = data2,
               mtry = mtr, 
               num.threads = 20, 
               probability = T, 
               importance = 'impurity',
               case.weights = data_w2,
               num.trees = ntrees,
               max.depth = mdepth)
  predic1 = predict(rf1, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  predic1 = predic1$predictions[,1] 
  predic2 = predict(rf2, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  predic2 = predic2$predictions[,1] 
  suit = data.frame(dt_test[,c("sourceTaxonName","targetTaxonName","interact")], predict_orig=predic1, predict_reduced=predic2)
  corrs <- ddply(suit, .(sourceTaxonName), split_corr)
  means <- ddply(suit, .(sourceTaxonName), split_mean)
  list(mean_cor=mean(corrs$COR),mean_change_prob=mean(means$mean_diff),num_removed = sel1)
}

#apply the function in a loop
set.seed(123)
output <- list()
for(i in 0:100){
  dt <- replicate(10,prdChange_TDTC_corrMOD(Int,Non,ins=2.5,out=1,abs=4.75,SD_foc,thresh=0.31,mtr=42,ntrees=400,mdepth=0,percPred = i))
  dt <- rbind(dt,i)
  output[[length(output) + 1]] <- dt  
  print(i)
}

saveRDS(output,"###/results/TDTC_preds_resultsChangeInt2Non_corr&mean_diff.rds")

#get the required data from the list together to plot
op <- data.frame(t(do.call("cbind", output)))
op <- data.frame(apply(op,2,unlist))
par(mfrow=c(2,2))
plot(op$i,op$mean_cor,xlab="% focal predator interact records remaining the same in training data", ylab = "correlation (pearson)", pch = 19, cex = 0.5, col=rgb(red=0.1, green=0.2, blue=0.2, alpha=0.3))
plot(op$i,op$mean_change_prob,xlab="% focal predator interact records remaining the same in training data", ylab = "mean_diff", pch = 19, cex = 0.5, col=rgb(red=0.1, green=0.2, blue=0.2, alpha=0.3))

