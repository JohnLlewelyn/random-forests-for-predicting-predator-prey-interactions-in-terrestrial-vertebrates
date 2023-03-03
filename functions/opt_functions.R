#opt functions

#prepare datasets
#Return also the weights (presences and absences sum to 1)
train_test_abs <- function(x,y,z,abs) {
  obs <- x
  obs$outside = NULL
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs <- y[sample(nrow(y), dim(obs)[1]*abs, replace = FALSE),] #unobs for training
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  test_obs = subset(data, sample == FALSE)
  test_obs = test_obs[test_obs$interact=="TRUE",]
  unobs2   = z
  unobs2   = unobs2[!(paste(unobs2$sourceTaxonName,unobs2$targetTaxonName)%in%paste(unobs$sourceTaxonName,unobs$targetTaxonName)),]
  unobs2   = unobs2[sample(nrow(unobs2),dim(test_obs)[1], replace = FALSE),]
  unobs2$interact <- as.factor(FALSE)
  unobs2$outside <- NULL
  test = rbind(test_obs,unobs2)
  train_w <- train$interact
  train_w <- ifelse(train_w=="TRUE",1,1/abs)
  return(list(train,test,train_w))
}

#identify columns with taxonomic information based on the presence of "Name"
tax_cols <- function (x) { which(grepl("Name",names(x)))}

#run RFs
optRF_thr <- function(x,y,z,abs,thresh,mtr,mns,ntrees,sampF,rpl,...) {
  ch <- train_test_abs(x,y,z,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  #scores = predic$predictions[,1]
  #lbls <- dt_test$interact
  #lbls <- ifelse(lbls=="TRUE",1,0)
  #all_auc <- auc(scores, lbls)
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  return(all_tss)
}

#get tss
#True skill statistic
tss <- function(pred_perf) {
  a <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"]
  d <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"]
  b <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"]
  c <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"]
  tss <- (a/(a+c))+(d/(b+d)) - 1
  return(tss)
}

#Matthews correlation coefficient
mcc <- function(pred_perf) {
  TP <- as.numeric(pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"])
  TN <- as.numeric(pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"])
  FP <- as.numeric(pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"])
  FN <- as.numeric(pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"])
  top <- (TP*TN)-(FP*FN)
  bottom <- sqrt((TP+FP)*(TP+FN)*(TN+FP)*(TN+FN))
  mcc <- top/bottom
  return(mcc)
}

#accuracy
acc <- function(pred_perf) {
  TP <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"]
  TN <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"]
  FP <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"]
  FN <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"]
acc <-  (TP+TN)/(TP+TN+FN+FP) 
return(acc)
}


#sensitivity/true positive rate/recall:	TP / (TP + FN)
sn <- function(pred_perf) {
  TP <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"]
  TN <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"]
  FP <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"]
  FN <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"]
  sn <-  TP/(TP+FN) 
  return(sn)
}

#specificity: SP	TN / (TN + FP)
sp <- function(pred_perf) {
  TP <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"]
  TN <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"]
  FP <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"]
  FN <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"]
  sp <-  TN/(TN+FP) 
  return(sp)
}

#run RFs and get multiple performance metrics
optRF_allM <- function(x,y,z,abs,thresh,mtr,mns,ntrees,sampF,rpl,...) {
  ch <- train_test_abs(x,y,z,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#run RFs and get multiple performance metrics with non-int selected based on size matching
optRF_allM_size <- function(x,y,z,abs,ins,out,thresh,mtr,mns,ntrees,sampF,rpl,...) {
  ch <- train_test_dfsw(x,y,z,abs,ins,out)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}


##train on GloBI+, but test on all SD
RF_SD <- function(x,y,dt_test, abs,thresh,mtr,mns,ntrees,sampF,rpl,...){
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs <- y[sample(nrow(y), dim(obs)[1]*abs, replace = FALSE),] 
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ .,  
               data = data,num.threads=20,probability = T,case.weights = data_w, mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl)
 
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}


#With selection of non-ints based on size
#With more sampled inside/outside range
train_test_dfsw <- function(x,y,z,ins,out,abs) { #ins and out = how many times more from inside versus outside suitable size range e.g., 10 and 1 or 1 and 10, abs = how many times more absent than present
  obs <- x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  #unobs$outside = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  test_obs = subset(data, sample == FALSE)
  test_obs = test_obs[test_obs$interact=="TRUE",]
  unobs2   = z
  unobs2   = unobs2[!(paste(unobs2$sourceTaxonName,unobs2$targetTaxonName)%in%paste(unobs$sourceTaxonName,unobs$targetTaxonName)),]
  unobs2   = unobs2[sample(nrow(unobs2),dim(test_obs)[1], replace = FALSE),]
  unobs2$interact <- as.factor(FALSE)
  #unobs2$outside = NULL
  test = rbind(test_obs,unobs2)
  test$outside = NULL
  train$outside <- ifelse(train$outside=="present",1,1/abs)
  train_w <- train$outside
  train$outside = NULL
  return(list(train,test,train_w))
}

#
#run RFs with size of non-int
optRF_SZ <- function(x,y,z,ins,out,abs,thresh,mtr,mns,ntrees,sampF,rpl,...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  #scores = predic$predictions[,1]
  #lbls <- dt_test$interact
  #lbls <- ifelse(lbls=="TRUE",1,0)
  #all_auc <- auc(scores, lbls)
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  return(all_tss)
}


##train on GloBI+, but test on all SD with int:out non-int size 
RF_SD_size <- function(x,y,dt_test, ins, out,abs,thresh,mtr,mns,ntrees,sampF,rpl,...){
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ .,  
               data = data,num.threads=20,probability = T,case.weights = data_w, mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl)
  
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

######function for connectance 0.2#######################################
#With more sampled inside/outside range
train_test_Sconn <- function(x,y,z,ins,out,abs,conn,wgt) { #conn = connectance - interaction/records in test data, wgt=how many times the combined weight of non-interactions compared to combined weight of interactions
  obs <- x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  #unobs$outside = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  test_obs = subset(data, sample == FALSE)
  test_obs = test_obs[test_obs$interact=="TRUE",]
  unobs2   = z
  unobs2   = unobs2[!(paste(unobs2$sourceTaxonName,unobs2$targetTaxonName)%in%paste(unobs$sourceTaxonName,unobs$targetTaxonName)),]
  unobs2   = unobs2[sample(nrow(unobs2),dim(test_obs)[1]*((1-conn)/conn), replace = FALSE),]
  unobs2$interact <- as.factor(FALSE)
  #unobs2$outside = NULL
  test = rbind(test_obs,unobs2)
  test$outside = NULL
  train$outside <- ifelse(train$outside=="present",1,(1/abs)*wgt)
  train_w <- train$outside
  train$outside = NULL
  return(list(train,test,train_w))
}

#run RFs with size of non-int, conn in testing set (int/(int+non)), and wgt for training set (combined weight of non/combined weight of int in training data)
optRF_Sconn <- function(x,y,z,ins,out,abs,conn, wgt, thresh,mtr,mns,ntrees,sampF,rpl,...) {
  ch <- train_test_Sconn(x,y,z,ins,out,abs,conn,wgt)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  return(all_tss)
}

#same as above but getting the other performance metrics
optRF_allM_Sconn <- function(x,y,z,ins,out,abs,conn, wgt, thresh,mtr,mns,ntrees,sampF,rpl,...) {
  ch <- train_test_Sconn(x,y,z,ins,out,abs,conn,wgt)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#same as above but train on global apply on Simpson Desert
RF_SD_Sconn <- function(x,y,dt_test, ins, out, abs, wgt,thresh,mtr,mns,ntrees,sampF,rpl,...){ #conn not required because not building testing set
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs*wgt)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ .,  
               data = data,num.threads=20,probability = T,case.weights = data_w, mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl)
  
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#same as above but adjust fewer parameters
RF_SD_Sconn_simp <- function(x,y,dt_test, ins, out, abs, wgt,thresh,mtr,ntrees,...){ #conn not required because not building testing set
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs*wgt)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ .,  
               data = data,num.threads=20,probability = T,case.weights = data_w, mtry = mtr, min.node.size = mns, 
               num.trees = ntrees, sample.fraction = sampF, replace = rpl)
  
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#run RFs with size of non-int, conn in testing set (int/(int+non)), and wgt for training set (combined weight of non/combined weight of int in training data)
#but adjust fewer parameters
optRF_Sconn_simp <- function(x,y,z,ins,out,abs,conn, wgt, thresh,mtr,ntrees,...) {
  ch <- train_test_Sconn(x,y,z,ins,out,abs,conn,wgt)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, 
               num.trees = ntrees) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  return(all_tss)
}

optRF_allM_Sconn_simp <- function(x,y,z,ins,out,abs,conn, wgt, thresh,mtr,ntrees,...) {
  ch <- train_test_Sconn(x,y,z,ins,out,abs,conn,wgt)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ .,  
               data = ch[[1]],num.threads=20,probability = T,case.weights = ch[[3]],mtry = mtr, 
               num.trees = ntrees) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}


#same as above but train on global apply on Simpson Desert
RF_SD_Sconn_simp <- function(x,y,dt_test, ins, out, abs, wgt,thresh,mtr,ntrees,...){ #conn not required because not building testing set
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs*wgt)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ .,  
               data = data,num.threads=20,probability = T,case.weights = data_w, mtry = mtr, 
               num.trees = ntrees)
  
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#just adjust depth, ins:out,abs,thresh,
optRF_depth <- function(x,y,z,ins,out,abs,thresh,mtr,mdepth,ntrees, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]],max.depth = mdepth,num.trees = ntrees) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  return(all_tss)
}
#same as above but with more performance metrics
optRF_allM_depth <- function(x,y,z,ins,out,abs,thresh,mtr,mdepth,ntrees, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]],max.depth = mdepth,num.trees = ntrees) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#and for applying above to Simpson Desert data
RF_SD_simp_depth <- function(x,y,dt_test, ins, out, abs,thresh,mtr,mdepth,ntrees,...){ #
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = data,mtry=mtr,num.threads=20,probability = T,case.weights = data_w, max.depth = mdepth, num.trees = ntrees) #
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  predic = predic$predictions[,1] > thresh
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}

#opt abs, thr, trees, depth######################
#train_test_abs <- function(x,y,z,abs) {
#just adjust depth,abs,thresh,trees
optRF_depth_NOsize <- function(x,y,z,abs,thresh,mtr,mdepth,ntrees, ...) {
  ch <- train_test_abs(x,y,z,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]],max.depth = mdepth,num.trees = ntrees) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  return(all_tss)
}

#same as above but for Simpson Desert
#and for applying above to Simpson Desert data
RF_SD_depth_NOsize <- function(x,y,dt_test, abs,thresh,mtr,mdepth,ntrees,...){ #
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs <- y[sample(nrow(y), dim(obs)[1]*abs, replace = FALSE),] 
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = data,mtry=mtr,num.threads=20,probability = T,case.weights = data_w, max.depth = mdepth, num.trees = ntrees) #
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  predic = predic$predictions[,1] > thresh
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_acc <- round(acc(pred_perf),3)
  all_mcc <- round(mcc(pred_perf),3)
  all_sn <- round(sn(pred_perf),3)
  all_sp <- round(sp(pred_perf),3)
  allM <- data.frame(cbind(c("tss","acc","mcc","sn","sp"),c(all_tss,all_acc,all_mcc,all_sn,all_sp)))
  return(list(allM))
}
