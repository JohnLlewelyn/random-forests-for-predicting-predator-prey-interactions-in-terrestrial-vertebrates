#Training dataset taxonomic coverage (TDTC) and performance - prey training dataset size (remove records)
#update file paths to run (search for lines with "###" to find where required)

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)

#get data
Int <- readRDS("###/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/allNon_sameCont.RDS")
SD_foc <- readRDS("###/allperms_cut2_20EVs.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0
SD_foc$source_aerial_mam <- 0

#get lists of Simpson Desert predators, prey (which has both)
prds <- unique(SD_foc$sourceTaxonName)
prey <- unique(SD_foc$targetTaxonName)

#remove cols
kp <- names(Non)[!(grepl("eig",names(Non)))]
Int <- Int[,names(Int)%in%c(kp,"interact","outside")]
Non <- Non[,names(Non)%in%c(kp,"interact","outside")]
SD_foc <- SD_foc[,names(SD_foc)%in%c(kp,"interact","outside")]

#the function for removing predators
prey_TDTC <- function(x,y,ins,out,abs,dt_test, thresh,mtr,percPrey,...){
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
  #split training data so SD_prey in one group
  SDprey <- data[data$targetTaxonName%in%prey,]
  Oprey <- data[!(data$targetTaxonName%in%prey),]
  #remove rows based on perc (remove from data and data_w)
  sel1 <- round(percPrey/100*length(SDprey$targetTaxonName))
  sel2 <- sample(1:length(SDprey$targetTaxonName),sel1)
  SDprey <- SDprey[sel2,]
  data <- rbind(SDprey,Oprey)
  data_w <- ifelse(data$interact=="TRUE",1,1/(table(data$interact)[2]/table(data$interact)[1]))
  dt_test = dt_test[,names(data)]
  dt_test$interact <- as.factor(dt_test$interact)
  data = data[,-tax_cols(data)]
  rf = ranger(formula = interact ~., 
              data = data,
              mtry = mtr, 
              num.threads = 20, 
              probability = T, 
              importance = 'impurity',
              case.weights = data_w)
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  lbls <- dt_test$interact
  lbls <- ifelse(lbls=="TRUE",1,0)
  all_auc <- auc(scores, lbls)
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  list(Tsk=c(all_tss),Ac=c(all_auc),num_removed = sel1)
}

#apply the function in a loop
output <- list()
for(i in 0:100){
  dt <- replicate(10,prey_TDTC(Int,Non,ins=1,out=1,abs=4,SD_foc,thresh=0.43,mtr=11,num.trees=800,max.depth=1000, percPrey = i))
  dt <- rbind(dt,i)
  output[[length(output) + 1]] <- dt   
}

saveRDS(output,"###/TDTC_prey_results.rds")

#get the required data from the list together
op <- data.frame(t(do.call("cbind", output)))
op <- data.frame(apply(op,2,unlist))

#plot it
plot(op$i,op$Tsk,xlab="% of records of focal prey remaining in training dataset", ylab = "TSS", pch = 19, cex = 0.5, col=rgb(red=0.1, green=0.2, blue=0.2, alpha=0.3))
