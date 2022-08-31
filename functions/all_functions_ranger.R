#functions for Ranger

#Fix column format, all as.character except for mass (numeric)
form_cols <- function (x) {
  x <- data.frame(sapply(x, as.character))
  x$sourceBodyMass.Value <- as.numeric(x$sourceBodyMass.Value)#scale(as.numeric(x$sourceBodyMass.Value))
  x$targetBodyMass.Value <- as.numeric(x$targetBodyMass.Value)#scale(as.numeric(x$targetBodyMass.Value))
  return(x)
}

#match observed interactions with unobserved (pseudo-absent) interactions, split into training and testing datasets by randomly splitting rows
#(allowing for training and testing pseudo-absenses to come from the same or different absences)
train_test_dfs <- function(x,y,z) {
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs <- y[sample(nrow(y), dim(obs)[1], replace = FALSE),] #unobs for training
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  if (identical(y,z)) {                                  #if same unobserved data used for training and testing, just split data according to sample
    test  = subset(data, sample == FALSE)
  } else {                                               #otherwise, keep the same observed interactions, and take unobserved interactions from another dataset
    test_obs = subset(data, sample == FALSE)
    test_obs = test_obs[test_obs$interact=="TRUE",]
    unobs2   = z
    unobs2   = unobs2[!(paste(unobs2$sourceTaxonName,unobs2$targetTaxonName)%in%paste(unobs$sourceTaxonName,unobs$targetTaxonName)),]
    unobs2   = unobs2[sample(nrow(unobs2),dim(test_obs)[1], replace = FALSE),]
    unobs2$interact <- as.factor(FALSE)
    test = rbind(test_obs,unobs2)
  }
  return(list(train,test))
}

#like  train_test_dfs but allow more pseudo-absences (abs = how many times more absent interactions). 
#Return also the weights (presences and absences sum to 1)
train_test_abs <- function(x,y,z,abs) {
  obs <- x
  obs$outside = NULL
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs <- y[sample(nrow(y), dim(obs)[1]*abs, replace = FALSE),] #unobs for training
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
  test = rbind(test_obs,unobs2)
  train_w <- train$interact
  train_w <- ifelse(train_w=="TRUE",1,1/abs)
  return(list(train,test,train_w))
}

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
  data <- rbind(obs,unobs[,names(obs)])
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

#With more sampled inside/outside range, twice as many absences in testing data
train_test_dfsw2 <- function(x,y,z,ins,out,abs) { #ins and out = how many times more from inside versus outside suitable size range e.g., 10 and 1 or 1 and 10, abs = how many times more absent than present
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
  unobs2   = unobs2[sample(nrow(unobs2),dim(test_obs)[1]*2, replace = FALSE),]
  unobs2$interact <- as.factor(FALSE)
  #unobs2$outside = NULL
  test = rbind(test_obs,unobs2)
  test$outside = NULL
  train$outside <- ifelse(train$outside=="present",1,1/abs)
  train_w <- train$outside
  train$outside = NULL
  return(list(train,test,train_w))
}


#With more sampled inside/outside range - split p10 predators from a particular continent out as the test dataset
train_test_dfs_cont <- function(x,y,ins,out,abs,prd_d,cont,num_pr) { #prd_d = preds and continent dist, cont is continent col name e.g., "sourceNorthAm"
  obs = x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),] #only include species with recorded predatory interactions
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  #unobs$outside = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs) #data is observed and unobserved, with ins:out and abs for unobserved
  preds = cont_sel(prd_d,cont,num_pr) #get the ## test predators
  train = subset(data, !(data$sourceTaxonName%in%preds)) #train using data without the focal predators
  test = subset(data, data$sourceTaxonName%in%preds)
  test <- test[test$interact=="TRUE",]   #take just the observed interactions of the test predators
  train$outside <- ifelse(train$outside=="present",1,1/abs)
  train_w <- train$outside
  train$outside = NULL
  test$outside = NULL
  pred_un <- y[y$sourceTaxonName%in%preds,] #now get all unobserved links for focal predator (so not ins:out abs)
  pred_un$outside = NULL
  pred_un$interact <- "FALSE"
  test <- rbind(test,pred_un)
  #each predator should have an equal number of observed and unobserved interactions
  test <- split(test,test$sourceTaxonName)
  #Sample the unobserved based on number of observed
  obc <- function(x) {x <- data.frame(x)
  x$obs=table(x$interact)[1]
  return(x)}
  test <- mapply(obc,test,SIMPLIFY = F)
  sel_unob <- function(x){
    x <- data.frame(x)
    x_ob <- x[x$interact=="TRUE",]
    x_unob <- x[x$interact=="FALSE",]
    x_unob <- x_unob[sample(nrow(x_unob),x$obs),]
    x <- rbind(x_ob,x_unob)
    return(x)
  }
  test <- mapply(sel_unob,test,SIMPLIFY = F)
  test <- do.call("rbind",test)
  return(list(train,test,train_w))
}

#With more sampled inside/outside range using all of GLobI for training, to test on another data set
train_test_allg <- function(x,y,ins,out,abs) { #ins and out = how many times more from inside versus outside suitable size range e.g., 10 and 1 or 1 and 10, abs = how many times more absent than present
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
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  return(list(data,data_w))
}



#With more sampled inside/outside range, and allowing more pseudo-absences than presences
train_test_dfswmore <- function(x,y,z,ins,out,num) { #ins and out = how many times more from inside versus outside suitable size range e.g., 10 and 1 or 1 and 10
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), dim(obs)[1]/(ins+out)*ins*num, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), dim(obs)[1]/(ins+out)*out*num, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$outside = NULL
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
  unobs2$outside = NULL
  test = rbind(test_obs,unobs2)
  return(list(train,test))
}


#match observed interactions with unobserved (pseudo-absent) interactions, split into training and testing datasets 
#by randomly allocating certain species to each data set
train_test_sp <- function(x,y,z) {  #x = observed interactions, y = unobserved interactions for training (all_unobs/size_unobs/fly_size_unob), z = unobserved interactions for testing (all_unob)
  obs <- x
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs <- y[sample(nrow(y), dim(obs)[1], replace = FALSE),]
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sp <- unique(data$targetTaxonName)
  rem <- sample(sp,0.25*length(sp))
  train = data[!(data$sourceTaxonName%in%rem), ]
  if (identical(y,z)) {
    test  = data[(data$sourceTaxonName%in%rem), ]
  } else {
    test_obs = data[(data$sourceTaxonName%in%rem)&data$interact=="TRUE", ]
    unobs2 = z
    unobs2$interact <- as.factor(FALSE)
    unobs2 = unobs2[!(paste(unobs2$sourceTaxonName,unobs2$targetTaxonName)%in%paste(train$sourceTaxonName,train$targetTaxonName)),]
    unobs2 = unobs2[sample(nrow(unobs2),dim(test_obs)[1], replace = FALSE),]
    test <- rbind(test_obs,unobs2)
  }
  return(list(train,test))
  #return(data)
}


#Prepare data so each predator has same number of absences as presences
train_test_sp_even <- function(x,y) {  #x = observed interactions, y = unobserved interactions 
  obs <- x
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  n <- as.vector(num_obs)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  split_unob <- split(y, y$sourceTaxonName)
  group_sizes <- vapply(split_unob, nrow, integer(1))
  sampled_obs <- mapply(sample, group_sizes, n)
  get_rows <- function(df, rows) df[rows, , drop = FALSE] 
  sampled_rows <- mapply(get_rows, split_unob, sampled_obs, SIMPLIFY = FALSE)
  unobs <- do.call(rbind, sampled_rows)
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  test = subset(data, sample == FALSE)
  return(list(train,test))
  
}

#Prepare data so absences:presences ratio the same for all predators - can't adjust inside:outside suitable prey size range for some predators, so skip this bit
train_test_even_abs <- function(x,y,abs) {  #x = observed interactions, y = unobserved interactions 
  obs <- x
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  num_df <- data.frame(sourceTaxonName=names(num_obs),obs = unname(num_obs))
  names(num_df) <- c("sourceTaxonName","rem","obs")
  num_df<- num_df[,c("sourceTaxonName","obs")]
  n <- as.vector(num_obs)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),] #only include observed predators as predators in the unobserved interactions
  y <- merge(y,num_df,by="sourceTaxonName", all.x=T) #add column indicating number of observed interactions
  split_unob <- split(y, y$sourceTaxonName)                 #split unobserved by predator species 
  group_sizes <- vapply(split_unob, nrow, integer(1)) #list of the number of rows of unobserved interactions for each predator
  sampled_unob <- mapply(sample, group_sizes, n*abs) #sample which rows to keep
  get_rows <- function(df, rows) df[rows, , drop = FALSE] 
  sampled_row <- mapply(get_rows, split_unob, sampled_unob, SIMPLIFY = FALSE) 
  unobs <- do.call(rbind, sampled_rows)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  unobs <- unobs[,names(obs)]
  data <- rbind(obs,unobs)
  #split by predator
  sp <- unique(data$targetTaxonName)
  rem <- sample(sp,0.25*length(sp))
  train = data[!(data$sourceTaxonName%in%rem), ]
  train_w <- ifelse(train$interact=="TRUE",1,1/abs)
  train$outside <- NULL
  test = data[(data$sourceTaxonName%in%rem), ]
  #make each predator in the test data have equal obs:unobs
  split_test <- split(test, test$sourceTaxonName)
  ev <- function(yy) {int = table(yy$interact)[which(names(table(yy$interact))=="TRUE")]
  yz = yy[yy$interact=="TRUE",]
  yw = yy[yy$interact=="FALSE",]
  yw = yw[sample(nrow(yw),int),]
  yz = rbind(yz,yw)
  return(yz)}
  split_test <- mapply(ev, split_test, SIMPLIFY = FALSE) 
  test <- do.call(rbind, split_test)
  return(list(train,test,train_w))
}

#Prepare data so each predator has same number of absences as presences - all GloBI for training before testing on SD
train_Gl_even <- function(x,y) {  #x = observed interactions, y = unobserved interactions 
  obs <- x
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  n <- as.vector(num_obs)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  split_unob <- split(y, y$sourceTaxonName)
  group_sizes <- vapply(split_unob, nrow, integer(1))
  sampled_obs <- mapply(sample, group_sizes, n)
  get_rows <- function(df, rows) df[rows, , drop = FALSE] 
  sampled_rows <- mapply(get_rows, split_unob, sampled_obs, SIMPLIFY = FALSE)
  unobs <- do.call(rbind, sampled_rows)
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  return(data)
}

#so can replicate training and testing on GloBI - randomly splitting data frame
rep_all <- function(x,y,z) {
  ch <- train_test_dfs(x,y,z)
  rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = ch[[1]],num.tree=500,mtry=5,num.threads=20) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred$predictions)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

#so can replicate training and testing on GloBI - randomly splitting data frame and including threshold
rep_all_thr <- function(x,y,z,thresh,mtr, ...) {
  ch <- train_test_dfs(x,y,z)
  rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],num.tree=500,mtry=mtr,num.threads=20,probability = T) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

rep_all_thrw <- function(x,y,z,ins,out,abs,thresh,mtr, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],num.tree=800,mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]]) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  #scores = predic$predictions[,1]
  #lbls <- dt_test$interact
  #lbls <- ifelse(lbls=="TRUE",1,0)
  #all_auc <- auc(scores, lbls)
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  #list(all_tss, all_auc)
  return(all_tss)
}
###same as above but don't specify mtr so can use default
rep_all_thrw_mtryD <- function(x,y,z,ins,out,abs,thresh, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],num.tree=500,num.threads=20,probability = T,case.weights = ch[[3]]) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  #scores = predic$predictions[,1]
  #lbls <- dt_test$interact
  #lbls <- ifelse(lbls=="TRUE",1,0)
  #all_auc <- auc(scores, lbls)
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  #list(all_tss, all_auc)
  return(all_tss)
}

rep_cont <- function(x,y,ins,out,abs,thresh,mtr,prd_d,cont,num_pr, ...) {
  ch <- train_test_dfs_cont(x,y,ins,out,abs, prd_d, cont,num_pr)  #x,y,z,ins,out,abs,prd_d,cont
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., data = ch[[1]],num.tree=500,mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]]) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

#so can specify trees usinf trees=, but redundant as ... includes num.trees argument
rep_all_thrw_trees <- function(x,y,z,ins,out,abs,thresh,mtr,trees, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],num.tree=trees,mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]]) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  #scores = predic$predictions[,1]
  #lbls <- dt_test$interact
  #lbls <- ifelse(lbls=="TRUE",1,0)
  #all_auc <- auc(scores, lbls)
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  #list(all_tss, all_auc)
  return(all_tss)
}

#so can replicate training and testing on GloBI - randomly splitting data frame by predatory species
rep_all_sp <- function(x,y,z,thresh,mtr, ...) {
  ch <- train_test_sp(x,y,z)
  rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],num.tree=500,mtry=mtr,num.threads=20,probability = TRUE) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

#so can replicate training and testing on GloBI - randomly splitting data frame
rep_all_ev <- function(x,y,thresh, mtr, ...) {
  ch <- train_test_sp_even(x,y)
  rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName,
               data = ch[[1]],num.tree=500,mtry=mtr,num.threads=20,probability = T) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

#so can replicate training on GloBI and testing on Simpson Desert
rep_SD <- function(x,y,z) {  # x = GloBI obs, y = GloBI unob, z = Simpson Desert (allperms)
  Gl <- train_Gl_SD(x,y)
  rf <- randomForest(interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = Gl,mtry=5, importance = TRUE) #so  
  pred_sd <- as.character(predict(rf, newdata=z[,-(which(names(z)=="interact"))], type = "class"))
  pred_sd_tab <- table(z[,which(names(z)=="interact")], pred_sd)
  all_tss <- round(tss(pred_sd_tab),3)
}


#GloBI data for training RF on before applying to other (e.g., Simpson Desert) datasets
train_Gl_SD <- function(x,y){
  obs <- x
  obs$interact <- as.factor(TRUE)
  unobs <- y[sample(nrow(y), dim(obs)[1], replace = FALSE),]
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  return(data)
}


#Use all observed GloBI + optimal number of unobserved links (abs) and optimal ratio inside and outside suitable size range (ins, out)
train_Gl_SDop <- function(x,y,ins,out,abs){
  ob <- x
  ob$interact <- as.factor(TRUE)
  unobs_ins <-y[y$outside=="FALSE",]
  unobs_ins <- unobs_ins[sample(nrow(unobs_ins), dim(ob)[1]*((ins/(ins+out))*abs), replace = FALSE),]
  unobs_out <-y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), dim(ob)[1]*((out/(ins+out))*abs), replace = FALSE),]
  unobs <- rbind(unobs_ins,unobs_out)
  unobs$interact <- as.factor(FALSE)
  unobs$outside = NULL
  data <- data.frame(rbind(ob,unobs))
  train_w <- c(ifelse(data$interact=="TRUE",1,1/abs))
  return(list(data,train_w))
}

#for making separate pseudo-absence sets to train and test
train_test_species <- function(x,y,z) {  #x = obs, y = size_unob/fly_size_unob, z = all_unob
  obs <- x
  obs$interact <- as.factor(TRUE)
  unobs <- y[sample(nrow(y), dim(obs)[1], replace = FALSE),] #matching number of unobserved from teh pseudo-absent set of interest
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sp <- unique(c(data$targetTaxonName,data$sourceTaxonName))
  rem <- sample(sp,0.3*length(sp))
  train <- data[!(data$sourceTaxonName%in%rem)&!(data$targetTaxonName%in%rem), ]
  test_ob  <- data[(data$sourceTaxonName%in%rem)&(data$targetTaxonName%in%rem)&data$interact=="TRUE", ]
  train_ab <- train[train$interact=="FALSE",]
  all_ab <- z[!(paste(z$sourceTaxonName,z$targetTaxonName)%in%paste(train_ab$sourceTaxonName,train_ab$targetTaxonName)),]
  all_ab$interact <- as.factor(FALSE)
  all_ab <- all_ab[sample(nrow(all_ab),dim(test_ob)[1],replace = FALSE),]
  test <- rbind(test_ob,all_ab)
  return(list(train,test))
  #return(data)
}


#True skill statistic
tss <- function(pred_perf) {
  a <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"]
  d <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"]
  b <- pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"]
  c <- pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"]
  tss <- (a/(a+c))+(d/(b+d)) - 1
  return(tss)
}

#True skill statistic to handle 1 column tables (e.g. when all predict the same)
tssF <- function(pred_perf) {
  a <- ifelse(length(pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"])==0,0,pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="TRUE"])
  d <- ifelse(length(pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"])==0,0,pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="FALSE"])
  b <- ifelse(length(pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"])==0,0,pred_perf[row.names(pred_perf)=="FALSE",colnames(pred_perf)=="TRUE"])
  c <- ifelse(length(pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"])==0,0,pred_perf[row.names(pred_perf)=="TRUE",colnames(pred_perf)=="FALSE"])
  tss <- (a/(a+c))+(d/(b+d)) - 1
  return(tss)
}

#observed and unobserved, without matching numbersies
train_test_allunob <- function(x,y) {
  obs <- x
  obs$interact <- as.factor(TRUE)
  unobs <- y
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  test  = subset(data, sample == FALSE)
  return(list(train,test))
  #return(data)
}


#get all permutations of two vectors
get.perms <- function(sdp,allsp){
  perms <- lapply(X=sdp, FUN=function(x,y) data.frame(cbind(x,y)), y = allsp)
  perms <-  do.call("rbind", perms)
  colnames(perms) <- c("sourceTaxonName","targetTaxonName")
  return(perms)
}

#plot rf predict results
plot_rf <- function(data, pred, main_T) {
  plot(log10(data$sourceBodyMass.Value),log10(data$targetBodyMass.Value),pch=19,col=as.factor(paste(pred,data$interact)),ylab= "log10(target mass)",xlab = "log10(source mass)", main = main_T,ylim=c(0,8))
}


#Function to change threshold in ranger and get TSS for predictions
rf_thresh <- function(dt_train,dt_test, thresh,mtr){
  rf = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, 
              data = dt_train,
              num.trees = 500, 
              mtry = mtr, 
              num.threads = 20, 
              probability = T, 
              importance = 'impurity')
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact")))])
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

#Function to change threshold in ranger and get TSS for predictions (and AUC), allowing weights for unobserved versus observed
rf_thresh_w <- function(x,y,ins,out,abs,dt_test, thresh,mtr,...){
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
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  #dt_test = dt_test[,!(tax_cols_sp(dt_test))]
  #try so there's an equal number observed vs unobserved in test data
  dt_ob <- dt_test[dt_test$interact=="TRUE",]
  dt_unob <- dt_test[dt_test$interact=="FALSE",]
  dt_unob <- dt_unob[sample(nrow(dt_unob),dim(dt_ob)[1]),]
  dt_test <- rbind(dt_ob,dt_unob)
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  #dt_test = dt_test[,!(names(dt_test)%in%dist)]
  rf = ranger(formula = interact ~., 
              data = data,
              num.trees = 500, 
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
  list(Tsk=c(all_tss),Ac=c(all_auc))
}

##SAME as above, train on GloBI+, but test on all SD
rf_thresh_all <- function(x,y,ins,out,abs,dt_test, thresh,mtr,...){
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
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
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
  list(Tsk=c(all_tss),Ac=c(all_auc))
}

#same as above no weights
rf_thresh_all_nw <- function(x,y,ins,out,abs,dt_test, thresh,mtr,...){
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
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
  rf = ranger(formula = interact ~., 
              data = data,
              mtry = mtr, 
              num.threads = 20, 
              probability = T, 
              importance = 'impurity',
              #case.weights = data_w
  )
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  lbls <- dt_test$interact
  lbls <- ifelse(lbls=="TRUE",1,0)
  all_auc <- auc(scores, lbls)
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  list(Tsk=c(all_tss),Ac=c(all_auc))
}


#Function to change threshold in ranger and get probability for each interaction, allowing weights for unobserved versus observed
rf_thresh_p <- function(x,y,ins,out,abs,dt_test,thresh,mtr){
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
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  #dt_test = dt_test[,!(tax_cols_sp(dt_test))]
  #try so there's an equal number observed vs unobserved in test data
  dt_ob <- dt_test[dt_test$interact=="TRUE",]
  dt_unob <- dt_test[dt_test$interact=="FALSE",]
  dt_unob <- dt_unob[sample(nrow(dt_unob),dim(dt_ob)[1]),]
  dt_test <- rbind(dt_ob,dt_unob)
  #dt_test = dt_test[,!(names(dt_test)%in%dist)]
  rf = ranger(formula = interact ~., 
              data = data,
              num.trees = 500, 
              mtry = mtr, 
              num.threads = 20, 
              probability = T, 
              importance = 'impurity',
              case.weights = data_w)
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","targetTaxonName","sourceTaxonName")))])
  predic = predic$predictions[,1]
  predic = data.frame(sourceTaxonName=dt_test$sourceTaxonName,targetTaxonName=dt_test$targetTaxonName,interact=dt_test$interact,prob=predic)
  return(list(predic))
}
#Class only taxanomic info = Function to change threshold in ranger and get TSS for predictions
rf_thresh_class <- function(dt_train,dt_test, thresh){
  rf = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName - sourceTaxonFamilyName - targetTaxonFamilyName - sourceTaxonOrderName - targetTaxonOrderName, 
              data = dt_train,
              num.trees = 500, 
              mtry = 5, 
              num.threads = 20, 
              probability = T)
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","targetTaxonName", "sourceTaxonName", "targetTaxonGenusName", "sourceTaxonGenusName", "sourceTaxonFamilyName", "targetTaxonFamilyName","sourceTaxonFamilyName", "targetTaxonFamilyName")))])
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}
#function to replicate rf_thresh, like in rep_all - dfs and sp
rep_thresh_dfs <- function(x,y,z,thresh) {
  ch <- train_test_dfs(x,y,z)
  rf_thresh(ch[[1]],ch[[2]],thresh)
  
}

rep_thresh_sp <- function(x,y,z,thresh) {
  ch <- train_test_sp(x,y,z)
  rf_thresh(ch[[1]],ch[[2]],thresh)
}
#and for SD
rep_thresh_SD <- function(x,y,z,thresh) {
  Gl <- train_Gl_SD(x,y)
  rf_thresh(Gl,z,thresh)
}

#Function to rerun RF #x and get threshhold for max TSS when splitting training and testing data by predatory species
best_thresh_sp <- function(obs,unob_train, unob_test,mtr){
  ch <- train_test_sp(obs,unob_train,unob_test)
  allG.ranger = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = ch[[1]], 
                       num.trees = 500, mtry = mtr, num.threads = 20, probability = TRUE)
  pred_prob = predict(allG.ranger, data = ch[[2]])$predictions
  #calculate tss for a given probability threshold - which probability threshold is best
  dt<-ch[[2]]
  iterations = 100
  variables = 2
  output <- matrix(ncol=variables,nrow=iterations)
  for(i in 1:iterations){
    pred_p <- factor(pred_prob[,1]>(i/100),levels = c("TRUE","FALSE"))
    predp_perf <- table(dt[,which(names(dt)=="interact")], pred_p)
    output[i,] <- c(round(tss(predp_perf),3),i/100)
    best_thresh <- output[which.max(output[,1]),2]
  }
  return(best_thresh)
}
#
#Function to rerun RF #x and get threshhold for max TSS when randomly splitting rows for training and testing data
best_thresh_dfs <- function(obs,unob_train, unob_test,mtr){
  ch <- train_test_dfs(obs,unob_train,unob_test)
  allG.ranger = ranger(formula = interact ~ .- targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName -sourceTaxonFamilyName -targetTaxonFamilyName -sourceTaxonOrderName -targetTaxonOrderName - sourceTaxonClassName -targetTaxonClassName, 
                       data = ch[[1]], num.trees = 500, mtry = mtr, num.threads = 20,importance = 'impurity',probability = T) 
  pred_prob = predict(allG.ranger, data = ch[[2]])$predictions
  #calculate tss for a given probability threshold
  dt<-ch[[2]]
  iterations = 100
  variables = 2
  output <- matrix(ncol=variables,nrow=iterations)
  for(i in 1:iterations){
    pred_p <- factor(pred_prob[,1]>(i/100),levels = c("TRUE","FALSE"))
    predp_perf <- table(dt[,which(names(dt)=="interact")], pred_p)
    output[i,] <- c(round(tss(predp_perf),3),i/100)
    best_thresh <- output[which.max(output[,1]),2]
  }
  return(best_thresh)
}


#Function to rerun RF #x and get threshhold for max TSS when matching species present vs absent data 50:50
best_thresh_ev <- function(obs,unob,mtr){
  ch <- train_test_sp_even(obs,unob)
  allG.ranger = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = ch[[1]], 
                       num.trees = 500, mtry = mtr, num.threads = 20, probability = TRUE)
  pred_prob = predict(allG.ranger, data = ch[[2]])$predictions
  #calculate tss for a given probability threshold - which probability threshold is best
  dt<-ch[[2]]
  iterations = 100
  variables = 2
  output <- matrix(ncol=variables,nrow=iterations)
  for(i in 1:iterations){
    pred_p <- factor(pred_prob[,1]>(i/100),levels = c("TRUE","FALSE"))
    predp_perf <- table(dt[,which(names(dt)=="interact")], pred_p)
    output[i,] <- c(round(tss(predp_perf),3),i/100)
    best_thresh <- output[which.max(output[,1]),2]
  }
  return(best_thresh)
}
#
##Function to specify threshold in ranger and get TSS for predictions for optimized model (inside:outside, present:absent)
rf_thresh_op <- function(dt_train,dt_test,train_w, thresh, mtr){
  rf = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, 
              data = dt_train,
              num.trees = 500, 
              mtry = mtr, 
              num.threads = 20, 
              probability = T,
              case.weights = train_w)
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact")))])
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}
#so can get optimised training data sets and apply fr_thresh_op in replicate
rep_thresh_op <- function(x,y,z,ins,out,abs,thresh) {
  ch = train_test_dfsw(x,y,z,ins,out,abs)
  rf_thresh_op(ch[[1]],ch[[2]],ch[[3]],thresh)
  
}

#Function to rerun RF #x and get threshhold for max TSS with inside:outside and present:absence optimised
best_thresh_op <- function(obs,unob_train, unob_test,ins,out,abs){
  ch <- train_test_dfsw(obs,unob_train,unob_test,ins,out,abs)
  allG.ranger = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = ch[[1]], 
                       num.trees = 500, mtry = 5, num.threads = 20, probability = TRUE, case.weights = ch[[3]])
  pred_prob = predict(allG.ranger, data = ch[[2]])$predictions
  #calculate tss for a given probability threshold
  dt<-ch[[2]]
  iterations = 100
  variables = 2
  output <- matrix(ncol=variables,nrow=iterations)
  for(i in 1:iterations){
    pred_p <- factor(pred_prob[,1]>(i/100),levels = c("TRUE","FALSE"))
    predp_perf <- table(dt[,which(names(dt)=="interact")], pred_p)
    output[i,] <- c(round(tss(predp_perf),3),i/100)
    best_thresh <- output[which.max(output[,1]),2]
  }
  return(best_thresh)
}

#
Gopt_SD <- function(obs, all_unobw, ins,out,abs, thresh){
  Gl <- train_Gl_SDop(obs, all_unobw, ins,out,abs) #prep data (15 in:5out, 1present:2absent)
  rf_thresh_op(Gl[[1]],allperms,Gl[[2]], thresh)} #rf

#function to get importance of different variable for multiple runs
imp_dfs <- function (x,y,z) {
  ch <- train_test_dfs(x,y,z) #
  rf <- ranger(formula = interact ~ ., 
               data = ch[[1]],
               num.trees = 500, 
               mtry = 5, 
               num.threads = 20,
               importance = 'impurity')
  return(importance(rf))}

imp_sp <- function (x,y,z) {
  ch <- train_test_sp(x,y,z) #
  rf <- ranger(formula = interact ~ ., 
               data = ch[[1]],
               num.trees = 500, 
               mtry = 5, 
               num.threads = 20,
               importance = 'impurity')
  return(importance(rf))}


#optimize mtry FOR DFS
#opt_mtry_dfs <- function(x,y,z) {   #obs_mr,all_unob_mr,all_unob_mr
#  ch <- train_test_dfs(x,y,z)#
#  train = ch[[1]]
#  test = ch[[2]]
#  output <- matrix(ncol=2,nrow=10)
#  for(i in 1:10) {
#    rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName -sourceTaxonFamilyName -targetTaxonFamilyName -sourceTaxonOrderName -targetTaxonOrderName - sourceTaxonClassName -targetTaxonClassName, data = train, num.trees = 500, mtry = i, num.threads = 20,importance = 'impurity')
#    pred <- predict(rf, data=test[,-(which(names(test)%in%c("interact","targetTaxonName","sourceTaxonName","targetTaxonGenusName", "sourceTaxonGenusName", "sourceTaxonFamilyName", "targetTaxonFamilyName", "sourceTaxonOrderName","targetTaxonOrderName","sourceTaxonClassName","targetTaxonClassName", data = train, num.trees = 500, mtry = i, num.threads = 20,importance = 'impurity')
#    pred_tab <- table(ch[[2]][,which(names(test)=="interact")], pred$predictions)
#    ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
#   output[i,] = c(round(tss(pred_tab),3),i) }
#   best_mtry <- output[which.max(output[,1]),2]
#   return(best_mtry)
# }

#optimize mtry FOR SP
opt_mtry_sp <- function(x,y,z){   #obs_mr,all_unob_mr,all_unob_mr
  ch <- train_test_sp(x,y,z)#
  train = ch[[1]]
  test = ch[[2]]
  output <- matrix(ncol=2,nrow=10)
  for(i in 1:10){
    rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = train, num.trees = 500, mtry = i, num.threads = 20,importance = 'impurity')
    pred <- predict(rf, data=test[,-(which(names(test)%in%c("interact")))])
    pred_tab <- table(ch[[2]][,which(names(test)=="interact")], pred$predictions)
    ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
    output[i,] = c(round(tss(pred_tab),3),i)}
  best_mtry <- output[which.max(output[,1]),2]
  return(best_mtry)
}

##optimize mtry FOR even
opt_mtry_ev <- function(x,y){   #obs_mr,all_unob_mr,all_unob_mr
  ch <- train_test_sp_even(x,y)#
  train = ch[[1]]
  test = ch[[2]]
  output <- matrix(ncol=2,nrow=10)
  for(i in 1:10){
    rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = train, num.trees = 500, mtry = i, num.threads = 20,importance = 'impurity')
    pred <- predict(rf, data=test[,-(which(names(test)%in%c("interact")))])
    pred_tab <- table(ch[[2]][,which(names(test)=="interact")], pred$predictions)
    ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
    output[i,] = c(round(tss(pred_tab),3),i)}
  best_mtry <- output[which.max(output[,1]),2]
  return(best_mtry)
}



##Function to rerun RF #x and get threshhold for max TSS randomly selecting absences for GloBI, and applying to SDs
best_thresh_SD <- function(obs,unob,mtr,SD,prep_train){
  ch <- prep_train(obs,unob)
  allG.ranger = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = ch, 
                       num.trees = 500, mtry = mtr, num.threads = 20, probability = TRUE)
  pred_prob = predict(allG.ranger, data = SD)$predictions
  #calculate tss for a given probability threshold - which probability threshold is best
  dt<-SD
  iterations = 100
  variables = 2
  output <- matrix(ncol=variables,nrow=iterations)
  for(i in 1:iterations){
    pred_p <- factor(pred_prob[,1]>(i/100),levels = c("TRUE","FALSE"))
    predp_perf <- table(dt[,which(names(dt)=="interact")], pred_p)
    output[i,] <- c(round(tss(predp_perf),3),i/100)
    best_thresh <- output[which.max(output[,1]),2]
  }
  return(best_thresh)
}


#train on GloBI and test on SD - even number of absences for each predator 
tt_even_SD <- function(globi_obs,globi_unob){
  obs <- globi_obs
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  n <- as.vector(num_obs)
  globi_unob <- globi_unob[globi_unob$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  split_unob <- split(globi_unob, globi_unob$sourceTaxonName)
  group_sizes <- vapply(split_unob, nrow, integer(1))
  sampled_obs <- mapply(sample, group_sizes, n)
  get_rows <- function(df, rows) df[rows, , drop = FALSE]
  sampled_rows <- mapply(get_rows, split_unob, sampled_obs, SIMPLIFY = FALSE)
  unobs <- do.call(rbind, sampled_rows)
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  return(data)
}


#Function to restrict unobserved interactions to species found on the same continent
same_cont <- function(x) {
  x <- x[x$sourceAfrica=="TRUE"&x$targetAfrica=="TRUE" | x$sourceEurope=="TRUE"&x$targetEurope=="TRUE" | x$sourceNorthAm=="TRUE"&x$targetNorthAm=="TRUE" | x$sourceSouthAm=="TRUE"&x$targetSouthAm=="TRUE" | x$sourceOceania=="TRUE"&x$targetOceania=="TRUE" | x$sourceAsia=="TRUE"&x$targetAsia=="TRUE" | x$sourceAntarctica=="TRUE"&x$targetAntarctica=="TRUE",]
  return(x)
}

#get the probabilities, with choosing mtry
rf_prob <- function(dt_train,dt_test,mtr){
  rf = ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, 
              data = dt_train,
              num.trees = 500, 
              mtry = mtr, 
              num.threads = 20, 
              probability = T)
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact")))])
  res <- data.frame(predic[[1]])
  res<- cbind(dt_test[,c("sourceTaxonName","targetTaxonName")],res)
  #colnames(res) <- c("yep","nope")
  return(res)}


#identify columns with taxonomic information based on teh presence of "Name"
tax_cols <- function (x) { which(grepl("Name",names(x)))}

#identify columns with taxonomic information based on the presence of "Name", 
tax_cols_sp <- function (x) { which(grepl("Name",names(x)[!(names(x)%in%c("sourceTaxonName","targetTaxonName"))]))}

#Function for just using the best represented target families
rep_fam_tst <- function(x,y,f){
  obs <- x[x$targetTaxonFamilyName%in%f,]
  obs$interact <- as.factor(TRUE)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  y <- y[y$targetTaxonFamilyName%in%f,]
  unobs <- y[sample(nrow(y), dim(obs)[1], replace = FALSE),] #unobs for training
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  sample = sample.split(data$interact, SplitRatio = .75)
  train = subset(data, sample == TRUE)
  test  = subset(data, sample == FALSE)
  rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, 
               data = train,
               num.trees = 500, 
               mtry = 6, 
               num.threads = 20,
               importance = 'impurity')
  pred <- predict(rf, data=test[,-(which(names(test)%in%c("interact")))])
  pred_tab <- table(test[,which(names(test)=="interact")], pred$predictions)
  ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
  return(round(tss(pred_tab),3))}

#function for selecting predator species from particular continents
cont_sel <- function(dt,cont,num_pr) {
  dd <- dt[dt[,cont]=="TRUE",]
  dd <- dd[sample(nrow(dd),num_pr),]
  return(dd$sourceTaxonName)
}

#function for selecting predator species from any continent
any_sel <- function(dt,num_pr) {
  dd=dt
  dd <- dd[sample(nrow(dd),num_pr),]
  return(dd$sourceTaxonName)
}

#With more sampled inside/outside range - split ## predators from any continent out as the test dataset
train_test_dfs_any <- function(x,y,ins,out,abs,prd_d,num_pr) { 
  obs = x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),] #only include species with recorded predatory interactions
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  #unobs$outside = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs) #data is observed and unobserved, with ins:out and abs for unobserved
  preds = any_sel(prd_d,num_pr) #get the ## test predators
  train = subset(data, !(data$sourceTaxonName%in%preds)) #train using data without the focal predators
  test = subset(data, data$sourceTaxonName%in%preds)
  test <- test[test$interact=="TRUE",]   #take just the observed interactions of the test predators
  train$outside <- ifelse(train$outside=="present",1,1/abs)
  train_w <- train$outside
  train$outside = NULL
  test$outside = NULL
  pred_un <- y[y$sourceTaxonName%in%preds,] #now get all unobserved links for focal predator (so not ins:out abs)
  pred_un$outside = NULL
  pred_un$interact <- "FALSE"
  test <- rbind(test,pred_un)
  #each predator should have an equal number of observed and unobserved interactions
  test <- split(test,test$sourceTaxonName)
  #Sample the unobserved based on number of observed
  obc <- function(x) {x <- data.frame(x)
  x$obs=table(x$interact)[1]
  return(x)}
  test <- mapply(obc,test,SIMPLIFY = F)
  sel_unob <- function(x){
    x <- data.frame(x)
    x_ob <- x[x$interact=="TRUE",]
    x_unob <- x[x$interact=="FALSE",]
    x_unob <- x_unob[sample(nrow(x_unob),x$obs),]
    x <- rbind(x_ob,x_unob)
    return(x)
  }
  test <- mapply(sel_unob,test,SIMPLIFY = F)
  test <- do.call("rbind",test)
  return(list(train,test,train_w))
}

#so can replicate any_sel
rep_any <- function(x,y,ins,out,abs,thresh,mtr,prd_d,num_pr, ...) {
  ch <- train_test_dfs_any(x,y,ins,out,abs, prd_d,num_pr)  #x,y,z,ins,out,abs,prd_d,cont
  rf <- ranger(formula = interact ~ . - targetTaxonName - sourceTaxonName - targetTaxonGenusName - sourceTaxonGenusName, data = ch[[1]],num.tree=500,mtry=mtr,num.threads=20,probability = T,case.weights = ch[[3]]) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  all_tss
}

#fix format of data so all columns, apart from names, are numeric or 1s and 0s
numfun <- function(x,y)  {   #where x is obs_mr and y is all_unobs_mro
  x$outside <- "observed"
  z <- rbind(x,y)
  z$sourceActivity.Nocturnal <- as.integer(z$sourceActivity.Nocturnal)
  z$sourceActivity.Diurnal <- as.integer(z$sourceActivity.Diurnal)
  z$sourceActivity.Crepuscular <- as.integer(z$sourceActivity.Crepuscular)
  z$sourceflight <- ifelse(z$sourceflight=="flies",1,0)
  z$sourceeat_plants <- as.integer(z$sourceeat_plants)
  z$sourceeat_verts <- as.integer(z$sourceeat_verts)
  z$sourceeat_invs <- as.integer(z$sourceeat_invs)
  z$targetActivity.Nocturnal <- as.integer(z$targetActivity.Nocturnal)
  z$targetActivity.Diurnal <- as.integer(z$targetActivity.Diurnal)
  z$targetActivity.Crepuscular <- as.integer(z$targetActivity.Crepuscular)
  z$targetflight <- ifelse(z$targetflight=="flies",1,0)
  z$targeteat_plants <- as.integer(z$targeteat_plants)
  z$targeteat_verts <- as.integer(z$targeteat_verts)
  z$targeteat_invs <- as.integer(z$targeteat_invs)
  z$source_marine_mam <-ifelse(z$sourceForStrat.Value=="M",1,0)
  z$source_ground_mam <-ifelse(z$sourceForStrat.Value=="G",1,0)
  z$source_scansorial_mam <-ifelse(z$sourceForStrat.Value=="S",1,0)
  z$source_arboreal_mam <-ifelse(z$sourceForStrat.Value=="Ar",1,0)
  z$source_aerial_mam <-ifelse(z$sourceForStrat.Value=="A",1,0)
  z$sourceForStrat.Value <- NULL
  z$target_marine_mam <-ifelse(z$targetForStrat.Value=="M",1,0)
  z$target_ground_mam <-ifelse(z$targetForStrat.Value=="G",1,0)
  z$target_scansorial_mam <-ifelse(z$targetForStrat.Value=="S",1,0)
  z$target_arboreal_mam <-ifelse(z$targetForStrat.Value=="Ar",1,0)
  z$target_aerial_mam <-ifelse(z$targetForStrat.Value=="A",1,0)
  z$targetForStrat.Value <- NULL
  z <- z[,-(which(names(z)%in%c("sourceTaxonGenusName","sourceTaxonFamilyName","sourceTaxonOrderName", "sourceTaxonClassName","targetTaxonGenusName","targetTaxonFamilyName","targetTaxonOrderName","targetTaxonClassName")))]
  obs = z[z$outside=="observed",]
  obs$outside = NULL
  unobs = z[z$outside!="observed",]
  return(list(obs,unobs))
}

#fix format of data so all columns, apart from names, are numeric or 1s and 0s for allperms
numfun_all <- function(z)  {   #where x is obs_mr and y is all_unobs_mro
  z$sourceActivity.Nocturnal <- as.integer(z$sourceActivity.Nocturnal)
  z$sourceActivity.Diurnal <- as.integer(z$sourceActivity.Diurnal)
  z$sourceActivity.Crepuscular <- as.integer(z$sourceActivity.Crepuscular)
  z$sourceflight <- ifelse(z$sourceflight=="flies",1,0)
  z$sourceeat_plants <- as.integer(z$sourceeat_plants)
  z$sourceeat_verts <- as.integer(z$sourceeat_verts)
  z$sourceeat_invs <- as.integer(z$sourceeat_invs)
  z$targetActivity.Nocturnal <- as.integer(z$targetActivity.Nocturnal)
  z$targetActivity.Diurnal <- as.integer(z$targetActivity.Diurnal)
  z$targetActivity.Crepuscular <- as.integer(z$targetActivity.Crepuscular)
  z$targetflight <- ifelse(z$targetflight=="flies",1,0)
  z$targeteat_plants <- as.integer(z$targeteat_plants)
  z$targeteat_verts <- as.integer(z$targeteat_verts)
  z$targeteat_invs <- as.integer(z$targeteat_invs)
  z$source_marine_mam <-ifelse(z$sourceForStrat.Value=="M",1,0)
  z$source_ground_mam <-ifelse(z$sourceForStrat.Value=="G",1,0)
  z$source_scansorial_mam <-ifelse(z$sourceForStrat.Value=="S",1,0)
  z$source_arboreal_mam <-ifelse(z$sourceForStrat.Value=="Ar",1,0)
  z$source_aerial_mam <-ifelse(z$sourceForStrat.Value=="A",1,0)
  z$sourceForStrat.Value <- NULL
  z$target_marine_mam <-ifelse(z$targetForStrat.Value=="M",1,0)
  z$target_ground_mam <-ifelse(z$targetForStrat.Value=="G",1,0)
  z$target_scansorial_mam <-ifelse(z$targetForStrat.Value=="S",1,0)
  z$target_arboreal_mam <-ifelse(z$targetForStrat.Value=="Ar",1,0)
  z$target_aerial_mam <-ifelse(z$targetForStrat.Value=="A",1,0)
  z$targetForStrat.Value <- NULL
  z <- z[,-(which(names(z)%in%c("sourceTaxonGenusName","sourceTaxonFamilyName","sourceTaxonOrderName", "sourceTaxonClassName","targetTaxonGenusName","targetTaxonFamilyName","targetTaxonOrderName","targetTaxonClassName")))]
  z <- z[,!(names(z)%in%dist)]
  return(z)
}


#Function to sample from inside and outside suitable size range
try_fun <- function(x,ins,out,abs){
  a <- x[x$outside=="TRUE",]
  b <- x[x$outside=="FALSE",]
  a <- a[sample(nrow(a),round(unique(x$obs)*abs*(out/(out+ins)))),]
  b <- b[sample(nrow(b),round(unique(x$obs)*abs*(ins/(out+ins)))),]
  xx <- rbind(a,b)
  return(as.data.frame(xx))
}

#function to sample unobserved rows, sampling from inside and outside is determined ratio if possible, otherwise samples at random across $outside
even_unob_samp <- function(x,ins,out,abs)  { 
  if(length(x$outside[x$outside=="TRUE"])<round((out/(out+ins))*abs*unique(x$obs))|length(x$outside[x$outside=="FALSE"])<round((ins/(ins+out))*abs*unique(x$obs))){
    y <- x[sample(nrow(x),unique(x$obs)*abs),]
  } else { 
    y <- try_fun(x,ins,out,abs)
  }
  return(y)
}


#Function to split GloBI data into training and testing, with obs_unobserved ratio same for all species, and controlling ins, out, abs, and allow for replication
even_rep <- function(x,y,ins,out,abs,mtr,thresh) {    #x <- obs_mr, y <- all_unob_mro
  obs <- x
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  num_df <- data.frame(sourceTaxonName=names(num_obs),obs = unname(num_obs))
  names(num_df) <- c("sourceTaxonName","rem","obs")
  num_df<- num_df[,c("sourceTaxonName","obs")]
  n <- as.vector(num_obs)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),] #only include observed predators as predators in the unobserved interactions
  y <- merge(y,num_df,by="sourceTaxonName", all.x=T) #add column indicating number of observed interactions
  split_unob <- split(y, y$sourceTaxonName)                 #split unobserved by predator species 
  group_sizes <- vapply(split_unob, nrow, integer(1)) #list of the number of rows of unobserved interactions for each predator
  sampled_unob <- mapply(even_unob_samp, split_unob,ins,out,abs,SIMPLIFY = FALSE) #sample which rows to keep
  unobs <- do.call(rbind, sampled_unob)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  unobs <- unobs[,names(obs)]
  data <- rbind(obs,unobs)
  #split by predator
  sp <- unique(data$targetTaxonName)
  rem <- sample(sp,0.25*length(sp))
  train = data[!(data$sourceTaxonName%in%rem), ]
  train_w <- ifelse(train$interact=="TRUE",1,1/abs)
  train$outside <- NULL
  test = data[(data$sourceTaxonName%in%rem), ]
  #make each predator in the test data have equal obs:unobs
  split_test <- split(test, test$sourceTaxonName)
  ev <- function(yy) {int = table(yy$interact)[which(names(table(yy$interact))=="TRUE")]
  yz = yy[yy$interact=="TRUE",]
  yw = yy[yy$interact=="FALSE",]
  yw = yw[sample(nrow(yw),int*abs),]
  yz = rbind(yz,yw)
  return(yz)}
  split_test <- mapply(ev, split_test, SIMPLIFY = FALSE) 
  test <- do.call(rbind, split_test)
  #remove all taxonomic information
  train <- train[,-tax_cols(train)]
  test <- test[,-tax_cols(test)]
  #run RF
  rf <- ranger(formula = interact ~ ., 
               data = train,
               num.trees = 500, 
               mtry = mtr, 
               num.threads = 20,
               case.weights = train_w,
               importance = 'impurity',
               probability = T,
               max.depth = 1000)
  #importance(rf) - mass of predator and prey is by far the most important
  pred <- predict(rf, data=test[,-(which(names(test)%in%c("interact")))])
  pred <- pred$predictions[,1] > thresh
  pred_tab <- table(test[,which(names(test)=="interact")], pred)
  ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
  return(round(tss(pred_tab),3)) 
}

#same as even_rep but use all GloBI for training, test on SD/allperms
even_repSD <- function(x,y,SD,ins,out,abs,mtr,thresh) {    #x <- int, y <- nono
  obs <- x
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  num_df <- data.frame(sourceTaxonName=names(num_obs),obs = unname(num_obs))
  names(num_df) <- c("sourceTaxonName","rem","obs")
  num_df<- num_df[,c("sourceTaxonName","obs")]
  n <- as.vector(num_obs)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),] #only include observed predators as predators in the unobserved interactions
  y <- merge(y,num_df,by="sourceTaxonName", all.x=T) #add column indicating number of observed interactions
  split_unob <- split(y, y$sourceTaxonName)                 #split unobserved by predator species 
  group_sizes <- vapply(split_unob, nrow, integer(1)) #list of the number of rows of unobserved interactions for each predator
  sampled_unob <- mapply(even_unob_samp, split_unob,ins,out,abs,SIMPLIFY = FALSE) #sample which rows to keep
  unobs <- do.call(rbind, sampled_unob)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  unobs <- unobs[,names(obs)]
  train <- rbind(obs,unobs)
  train_w <- ifelse(train$interact=="TRUE",1,1/abs)
  train$outside <- NULL
  test = SD
  #make each predator in the test data have equal obs:unobs
  split_test <- split(test, test$sourceTaxonName)
  ev <- function(yy) {int = table(yy$interact)[which(names(table(yy$interact))=="TRUE")]
  yz = yy[yy$interact=="TRUE",]
  yw = yy[yy$interact=="FALSE",]
  yw = yw[sample(nrow(yw),int),]
  yz = rbind(yz,yw)
  return(yz)}
  split_test <- mapply(ev, split_test, SIMPLIFY = FALSE) 
  test <- do.call(rbind, split_test)
  #remove all taxonomic information
  train <- train[,-tax_cols(train)]
  test <- test[,-tax_cols(test)]
  #run RF
  rf <- ranger(formula = interact ~ ., 
               data = train,
               num.trees = 500, 
               mtry = mtr, 
               num.threads = 20,
               case.weights = train_w,
               importance = 'impurity',
               probability = T,
               max.depth = 1000)
  #importance(rf) - mass of predator and prey is by far the most important
  pred <- predict(rf, data=test[,-(which(names(test)%in%c("interact")))])
  pred <- pred$predictions[,1] > thresh
  pred_tab <- table(test[,which(names(test)=="interact")], pred)
  ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
  return(round(tss(pred_tab),3)) 
}


#function to get probabilities with even obs:unobs ratios between predators
even_repSDp <- function(x,y,SD,ins,out,abs,mtr,thresh) {    #x <- int, y <- nono
  obs <- x
  obs$interact <- as.factor(TRUE)
  num_obs <- table(obs$sourceTaxonName) #how many observation per predator
  num_df <- data.frame(sourceTaxonName=names(num_obs),obs = unname(num_obs))
  names(num_df) <- c("sourceTaxonName","rem","obs")
  num_df<- num_df[,c("sourceTaxonName","obs")]
  n <- as.vector(num_obs)
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),] #only include observed predators as predators in the unobserved interactions
  y <- merge(y,num_df,by="sourceTaxonName", all.x=T) #add column indicating number of observed interactions
  split_unob <- split(y, y$sourceTaxonName)                 #split unobserved by predator species 
  group_sizes <- vapply(split_unob, nrow, integer(1)) #list of the number of rows of unobserved interactions for each predator
  sampled_unob <- mapply(even_unob_samp, split_unob,ins,out,abs,SIMPLIFY = FALSE) #sample which rows to keep
  unobs <- do.call(rbind, sampled_unob)
  unobs$interact <- as.factor(FALSE)
  unobs$outside <- NULL
  unobs <- unobs[,names(obs)]
  train <- rbind(obs,unobs)
  train_w <- ifelse(train$interact=="TRUE",1,1/abs)
  train$outside <- NULL
  test = SD
  #make each predator in the test data have equal obs:unobs
  split_test <- split(test, test$sourceTaxonName)
  ev <- function(yy) {int = table(yy$interact)[which(names(table(yy$interact))=="TRUE")]
  yz = yy[yy$interact=="TRUE",]
  yw = yy[yy$interact=="FALSE",]
  yw = yw[sample(nrow(yw),int),]
  yz = rbind(yz,yw)
  return(yz)}
  split_test <- mapply(ev, split_test, SIMPLIFY = FALSE) 
  test <- do.call(rbind, split_test)
  #remove all taxonomic information
  train <- train[,-tax_cols(train)]
  test <- test[,-tax_cols(test)]
  #run RF
  rf <- ranger(formula = interact ~ ., 
               data = train,
               num.trees = 500, 
               mtry = mtr, 
               num.threads = 20,
               case.weights = train_w,
               importance = 'impurity',
               probability = T,
               max.depth = 1000)
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","targetTaxonName","sourceTaxonName")))])
  predic = predic$predictions[,1]
  predic = data.frame(sourceTaxonName=dt_test$sourceTaxonName,targetTaxonName=dt_test$targetTaxonName,interact=dt_test$interact,prob=predic)
  return(list(predic)) }

#function to best ins and outs
repfun <-function(obs, unob,ins,outs,abs) {
  ch <- train_test_dfsw(obs,unob,unob,ins,outs,abs) #
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., 
               data = ch[[1]],
               num.trees = 500, 
               mtry = 9, 
               num.threads = 20,
               case.weights = ch[[3]],
               probability = T)
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact")))])
  pred <- pred$predictions[,1]>0.49
  pred_tab <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  ((pred_tab[1,2] + pred_tab[2,1])/sum(pred_tab))*100  #OOB estimate of  error rate
  #output[ins,] <- c(round(tss(pred_tab),3),ins)}
  round(tss(pred_tab),3)}

#AUC function
auc <- function( scores, lbls ){
  stopifnot( length(scores) == length(lbls) )
  jp <- which( lbls > 0 ); np <- length( jp )
  jn <- which( lbls <= 0); nn <- length( jn )
  s0 <- sum( rank(scores)[jp] )
  (s0 - np*(np+1) / 2) / (np*nn)
}   

#add phylogenetic information from RDS files to GloBI (consM and rawM)
addPhy <- function(phy,ob,un){
  vecs <- phy[[1]]
  dist <- phy[[2]]
  vecs <- vecs[,1:3]
  tvecs <- svecs <- vecs
  names(tvecs) <- c("targetEV1","targetEV2","targetEV3")
  names(svecs) <- c("sourceEV1","sourceEV2","sourceEV3")
  tvecs$targetTaxonName<-gsub("_", " ",row.names(tvecs))
  svecs$sourceTaxonName<-gsub("_", " ",row.names(svecs))
  dist <- melt(dist)
  names(dist) <- c("targetTaxonName","sourceTaxonName","phyDist")
  dist$targetTaxonName<- gsub("_", " ",dist$targetTaxonName)
  dist$sourceTaxonName<- gsub("_", " ",dist$sourceTaxonName)
  
  #now merge the data, eigenvectors first
  ob <- merge(ob,tvecs, by="targetTaxonName")
  ob <- merge(ob,svecs, by="sourceTaxonName")
  un <- merge(un,tvecs, by="targetTaxonName")
  un <- merge(un,svecs, by="sourceTaxonName")
  #now distances
  ob$x = paste(ob$targetTaxonName,ob$sourceTaxonName)
  dist$x = paste(dist$targetTaxonName,dist$sourceTaxonName)
  un$x = paste(un$targetTaxonName,un$sourceTaxonName)
  dist$targetTaxonName<-dist$sourceTaxonName<-NULL
  ob <- merge(ob, dist, by = "x")
  un <- merge(un, dist, by = "x")
  ob$x<-dist$x<-un$x<-NULL
  return(list(ob,un))
}

#add phylogenetic information from RDS files to SD (consM and rawM)
addPhySD <- function(phy,SD){
  vecs <- phy[[1]]
  dist <- phy[[2]]
  vecs <- vecs[,1:3]
  tvecs <- svecs <- vecs
  names(tvecs) <- c("targetEV1","targetEV2","targetEV3")
  names(svecs) <- c("sourceEV1","sourceEV2","sourceEV3")
  tvecs$targetTaxonName<-gsub("_", " ",row.names(tvecs))
  svecs$sourceTaxonName<-gsub("_", " ",row.names(svecs))
  dist <- melt(dist)
  names(dist) <- c("targetTaxonName","sourceTaxonName","phyDist")
  dist$targetTaxonName<- gsub("_", " ",dist$targetTaxonName)
  dist$sourceTaxonName<- gsub("_", " ",dist$sourceTaxonName)
  
  #now merge the data, eigenvectors first
  SD <- merge(SD,tvecs, by="targetTaxonName")
  SD <- merge(SD,svecs, by="sourceTaxonName")
  
  #now distances
  SD$x = paste(SD$targetTaxonName,SD$sourceTaxonName)
  dist$x = paste(dist$targetTaxonName,dist$sourceTaxonName)
  dist$targetTaxonName<-dist$sourceTaxonName<-NULL
  SD <- merge(SD, dist, by = "x")
  SD$x<-dist$x<-NULL
  return(SD)
}

rep_all_thrnw <- function(x,y,z,ins,out,abs,thresh,mtr, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., # - targetTaxonGenusName - sourceTaxonGenusName, 
               data = ch[[1]],num.tree=500,mtry=mtr,num.threads=20,probability = T) #
  pred <- predict(rf, data=ch[[2]][,-(which(names(ch[[2]])%in%c("interact","sourceTaxonName","targetTaxonName")))])
  #scores = predic$predictions[,1]
  #lbls <- dt_test$interact
  #lbls <- ifelse(lbls=="TRUE",1,0)
  #all_auc <- auc(scores, lbls)
  pred = pred$predictions[,1] > thresh
  pred_perf <- table(ch[[2]][,which(names(ch[[2]])=="interact")], pred)
  all_tss <- round(tss(pred_perf),3)
  #list(all_tss, all_auc)
  return(all_tss)
}

#Function to change threshold in ranger and get TSS for predictions (and AUC), allowing weights for unobserved versus observed
rf_thresh_nw <- function(x,y,ins,out,abs,dt_test, thresh,mtr,...){
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
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  #dt_test = dt_test[,!(tax_cols_sp(dt_test))]
  #try so there's an equal number observed vs unobserved in test data
  dt_ob <- dt_test[dt_test$interact=="TRUE",]
  dt_unob <- dt_test[dt_test$interact=="FALSE",]
  dt_unob <- dt_unob[sample(nrow(dt_unob),dim(dt_ob)[1]),]
  dt_test <- rbind(dt_ob,dt_unob)
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  #dt_test = dt_test[,!(names(dt_test)%in%dist)]
  rf = ranger(formula = interact ~., 
              data = data,
              num.trees = 500, 
              mtry = mtr, 
              num.threads = 20, 
              probability = T, 
              importance = 'impurity',
              #case.weights = data_w
  )
  predic = predict(rf, data=dt_test[,-(which(names(dt_test)%in%c("interact","sourceTaxonName","targetTaxonName")))])
  scores = predic$predictions[,1]
  lbls <- dt_test$interact
  lbls <- ifelse(lbls=="TRUE",1,0)
  all_auc <- auc(scores, lbls)
  predic = predic$predictions[,1] > thresh
  #res <- data.frame(cbind(as.logical(dt_test$interact),predic))
  pred_perf <- table(dt_test$interact, predic)
  all_tss <- round(tss(pred_perf),3)
  list(Tsk=c(all_tss),Ac=c(all_auc))
}

#Function to change threshold in ranger and get TSS for predictions (and AUC), allowing weights for unobserved versus observed, and specifically including set number of unobserved from focal predators
rf_thresh_w_SDnon <- function(x,y,ins,out,abs,dt_test, thresh,mtr,prds,...){
  obs <- x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  SD_num_ob <- dim(obs[obs$sourceTaxonName%in%prds,])[1]#how many observed interaction for focal predator from inside
  unobs_inSD <- unobs_in[unobs_in$sourceTaxonName%in%prds,] #get all the NonInteractions for focal predators
  unobs_inSD <- unobs_inSD[sample(nrow(unobs_inSD),(SD_num_ob/(ins+out)*ins)*abs, replace=FALSE),]#select the NonInteractions for SD preds to keep
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs-dim( unobs_inSD)[1], replace = FALSE),] #other unobs for adding
  unobs_in <- rbind(unobs_in,unobs_inSD) #combine the SD and other 'inside suitable prey size' data
  unobs_out <- y[y$outside=="TRUE",]
  unobs_outSD <- unobs_out[unobs_out$sourceTaxonName%in%prds,] #get all the NonInteractions for focal predators
  unobs_outSD <- unobs_outSD[sample(nrow(unobs_outSD),(SD_num_ob/(ins+out)*out)*abs, replace=FALSE),]#select the NonInteractions for SD preds to keep
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs-dim( unobs_outSD)[1], replace = FALSE),]
  unobs_out <- rbind(unobs_out,unobs_outSD)
  unobs <- rbind(unobs_in,unobs_out)
  unobs$obs = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  #dt_test = dt_test[,!(tax_cols_sp(dt_test))]
  #try so there's an equal number observed vs unobserved in test data
  dt_ob <- dt_test[dt_test$interact=="TRUE",]
  dt_unob <- dt_test[dt_test$interact=="FALSE",]
  dt_unob <- dt_unob[sample(nrow(dt_unob),dim(dt_ob)[1]),]
  dt_test <- rbind(dt_ob,dt_unob)
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  #dt_test = dt_test[,!(names(dt_test)%in%dist)]
  rf = ranger(formula = interact ~., 
              data = data,
              num.trees = 500, 
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
  list(Tsk=c(all_tss),Ac=c(all_auc))
}

##SAME as above, train on GloBI+, but test on all SD
rf_thresh_all_SDnon <- function(x,y,ins,out,abs,dt_test, thresh,mtr,prds,...){
  obs <- x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  SD_num_ob <- dim(obs[obs$sourceTaxonName%in%prds,])[1]#how many observed interaction for focal predator from inside
  unobs_inSD <- unobs_in[unobs_in$sourceTaxonName%in%prds,] #get all the NonInteractions for focal predators
  unobs_inSD <- unobs_inSD[sample(nrow(unobs_inSD),(SD_num_ob/(ins+out)*ins)*abs, replace=FALSE),]#select the NonInteractions for SD preds to keep
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs-dim( unobs_inSD)[1], replace = FALSE),] #other unobs for adding
  unobs_in <- rbind(unobs_in,unobs_inSD) #combine the SD and other 'inside suitable prey size' data
  unobs_out <- y[y$outside=="TRUE",]
  unobs_outSD <- unobs_out[unobs_out$sourceTaxonName%in%prds,] #get all the NonInteractions for focal predators
  unobs_outSD <- unobs_outSD[sample(nrow(unobs_outSD),(SD_num_ob/(ins+out)*out)*abs, replace=FALSE),]#select the NonInteractions for SD preds to keep
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs-dim( unobs_outSD)[1], replace = FALSE),]
  unobs_out <- rbind(unobs_out,unobs_outSD)
  unobs <- rbind(unobs_in,unobs_out)
  unobs$obs = NULL
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  data_w <- ifelse(data$interact=="TRUE",1,1/abs)
  data$outside = NULL
  data = data[,-tax_cols(data)]  #ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  dt_test = dt_test[,names(data)]
  #dt_test <- dt_test[,-tax_cols(dt_test)] #can include species names or remove them, doesn't make a difference
  dt_test$interact <- as.factor(dt_test$interact)
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
  list(Tsk=c(all_tss),Ac=c(all_auc))
}

#replace original 3 eigenvectors (PEMs) with 20 evs
replaceEV <- function(x,ev,i){
  x <- x[,!(grepl("EV", names(x)))]
  x <- x[,!(grepl("phyDist", names(x)))]
  ev <- ev[[i]][[1]]
  ev$TaxonName <- row.names(ev)
  tvecs <- svecs <- ev
  names(tvecs) <- paste("target", names(tvecs), sep="")
  names(svecs) <- paste("source", names(svecs), sep = "")
  tvecs$targetTaxonName<-gsub("_", " ",row.names(tvecs))
  svecs$sourceTaxonName<-gsub("_", " ",row.names(svecs))
  tbb <- tvecs[tvecs$targetTaxonName=="Bubalus depressicornis",]
  sbb <- svecs[svecs$sourceTaxonName=="Bubalus depressicornis",]
  tbb$targetTaxonName="Bubalus bubalis"
  sbb$sourceTaxonName="Bubalus bubalis"
  tvecs <- rbind(tvecs,tbb)
  svecs <- rbind(svecs, sbb)
  #now merge the data, eigenvectors first
  x <- merge(x,tvecs, by="targetTaxonName")
  x <- merge(x,svecs, by="sourceTaxonName")
  return(x)
}

#Get importance of variables so can avoid over-fitting
imp_var <- function (x,y,z,ins,out,abs){
ch <- train_test_dfsw(x,y,z,ins,out,abs)
ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
rf <- ranger(formula = interact ~ ., 
             data = ch[[1]],
             num.trees = 600, 
             mtry = 1, 
             num.threads = 20,
             case.weights = ch[[3]],
             probability = T,
             splitrule = "gini", #not much difference, but "gini' appears to do a bit better than "hellinger"
             #max.depth = 25,  #can optimise this variable
             importance = 'permutation')
imp <- data.frame(importance(rf)) 
return(imp)}

all_imp <- function(x,y,z,ins,out,abs,thresh,mtr, ...) {
  ch <- train_test_dfsw(x,y,z,ins,out,abs)
  ch[[1]] <- ch[[1]][,-tax_cols(ch[[1]])]
  ch[[2]] <- ch[[2]][,-tax_cols(ch[[2]])]
  rf <- ranger(formula = interact ~ ., 
               data = ch[[1]],
               num.tree=600,
               mtry=mtr,
               num.threads=20,
               probability = T,
               case.weights = ch[[3]],
               importance = 'permutation') #
  imp <- data.frame(importance(rf)) 
  return(imp)
}

