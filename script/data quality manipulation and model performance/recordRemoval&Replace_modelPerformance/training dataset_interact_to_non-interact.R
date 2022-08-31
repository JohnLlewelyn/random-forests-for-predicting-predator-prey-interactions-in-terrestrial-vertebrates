#Training dataset and false non-interactions (switch interaction to non-interactions) in training data
#update file paths to run (search for lines with "###" to find where required)

setwd("~/###")
source("###/all_functions_ranger.R")

#libraries
library(ranger)

#get the data
Int <- readRDS("###/GloBIplus_Int20EVs.RDS")
Non <- readRDS("###/allNon_sameCont.RDS")
SD_foc <- readRDS("###/allperms_cut2_20EVs.RDS")

#add source_aerial_mam column to Int (because it's in the target and therefore potentially the noninteraction source column)
Int$source_aerial_mam <- 0
SD_foc$source_aerial_mam <- 0

#select columns to keep
nms <- names(Non)[!(grepl("eig",names(Non)))]
Int <- Int[,names(Int)%in%c(nms,"interact","outside")]
Non <- Non[,names(Non)%in%c(nms,"interact","outside")]
SD_foc <- SD_foc[,names(SD_foc)%in%c(nms,"interact","outside")]

#the function
rf_TDS <- function(x,y,ins,out,abs,dt_test, thresh,mtr,percI,...){
  obs <- x
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
  dt_test = dt_test[,names(data)]
  #switch interact rows to noninteraction rows based on perc (remove from data and data_w)
  Iprd <- data[data$interact==TRUE,]
  Nprd <- data[data$interact==FALSE,]
  sel1 <- round(percI/100*length(Iprd$sourceTaxonName))
  sel2 <- sample(1:length(Iprd$sourceTaxonName),sel1)
  Iprd$interact[sel2] <- FALSE
  data <- rbind(Iprd,Nprd)
  data_w <- ifelse(data$interact=="TRUE",1,1/(table(data$interact)[2]/table(data$interact)[1]))
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  dt_test$interact <- as.factor(dt_test$interact)
  dt_test <- dt_test[,names(data)]
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
  all_tss <- round(tssF(pred_perf),3)
  round(tss(pred_perf),3)
  round(tssF(pred_perf),3)
  list(Tsk=c(all_tss),Ac=c(all_auc),num_removed = sel1)
}


#put it in a for loop
output <- list()
for(i in 1:100){
  dt <- replicate(10,rf_TDS(Int,Non,ins=1,out=1,abs=4,SD_foc,thresh=0.43,mtr=11,num.trees=800,max.depth=1000, perc = i))
  dt <- rbind(dt,i)
  output[[length(output) + 1]] <- dt  
  print(i)
}
saveRDS(output,"###/falseNeg_anySpecies_results.rds")

#get the required data from the list together
op <- data.frame(t(do.call("cbind", output)))
op <- data.frame(apply(op,2,unlist))
op$i <- 100-op$i

plot(op$i,op$Tsk,xlab="% of interactions training dataset remaning as interactions", ylab = "TSS", pch = 19, cex = 0.5, col=rgb(red=0.1, green=0.2, blue=0.2, alpha=0.3))
