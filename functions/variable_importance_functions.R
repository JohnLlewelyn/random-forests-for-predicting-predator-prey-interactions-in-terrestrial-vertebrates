#functions for identifying important variables

#function to create one dataset from Ints and Non, specifying ratio of non-interactions from within preferred size range, and ratio of non-interactions to interactions
RF_allInt <- function(x,y,ins,out,abs) { 
  obs <- x
  obs$interact <- as.factor(TRUE)
  obs$outside = "present"
  y <- y[y$sourceTaxonName%in%unique(obs$sourceTaxonName),]
  unobs_in <- y[y$outside=="FALSE",]
  unobs_in <- unobs_in[sample(nrow(unobs_in), (dim(obs)[1]/(ins+out)*ins)*abs, replace = FALSE),] #unobs for training
  unobs_out <- y[y$outside=="TRUE",]
  unobs_out <- unobs_out[sample(nrow(unobs_out), (dim(obs)[1]/(ins+out)*out)*abs, replace = FALSE),]
  unobs <- rbind(unobs_in,unobs_out)
  unobs$interact <- as.factor(FALSE)
  data <- rbind(obs,unobs)
  data$outside <- ifelse(data$outside=="present",1,1/abs)
  data_w <- data$outside
  data$outside <- NULL
  data = data[,-tax_cols(data)]  
  return(list(data,data_w))
}

#identify columns with taxonomic information based on the presence of "Name"
tax_cols <- function (x) { which(grepl("Name",names(x)))}

#function to identify the least important variable to an RF, taking the mean importance score across source and target species
least_imp <- function(dat, wgt, mtryv, ntreesv, depthv) {
  rf <- ranger(formula = interact ~ ., 
               data = dat,
               num.trees = ntreesv, 
               mtry = mtryv, 
               num.threads = 20,
               case.weights = wgt,
               probability = T,
               splitrule = "gini", #not much difference, but "gini' appears to do a bit better than "hellinger"
               max.depth = depthv,
               importance = 'permutation')
  imp <- data.frame(importance(rf)) 
  #remove "source" and "target" so can aggregate importance
  imp$var<- gsub("source","", gsub("target", "", rownames(imp)))
  #aggregate
  agg <- aggregate(imp$importance.rf.,by=list(imp$var),FUN=mean)
  names(agg) <- c("var","importance")
  agg <- agg[order(agg$importance,decreasing = TRUE),]
  return(agg$var[which.min(agg$importance)])}

#least_imp but using default hyper parameters
least_imp_default <- function(dat, wgt) {
  rf <- ranger(formula = interact ~ ., 
               data = dat,
               num.threads = 20,
               case.weights = wgt,
               probability = T,
               splitrule = "gini", #not much difference, but "gini' appears to do a bit better than "hellinger"
               importance = 'permutation')
  imp <- data.frame(importance(rf)) 
  #remove "source" and "target" so can aggregate importance
  imp$var<- gsub("source","", gsub("target", "", rownames(imp)))
  #aggregate
  agg <- aggregate(imp$importance.rf.,by=list(imp$var),FUN=mean)
  names(agg) <- c("var","importance")
  agg <- agg[order(agg$importance,decreasing = TRUE),]
  return(agg$var[which.min(agg$importance)])}

#least_imp but using default hyper parameters and only considering phylo eigenvectors
least_imp_default_eig <- function(dat, wgt) {
  rf <- ranger(formula = interact ~ ., 
               data = dat,
               num.threads = 20,
               case.weights = wgt,
               probability = T,
               splitrule = "gini", #not much difference, but "gini' appears to do a bit better than "hellinger"
               importance = 'permutation')
  imp <- data.frame(importance(rf)) 
  #remove "source" and "target" so can aggregate importance
  imp$var<- gsub("source","", gsub("target", "", rownames(imp)))
  #aggregate
  agg <- aggregate(imp$importance.rf.,by=list(imp$var),FUN=mean)
  names(agg) <- c("var","importance")
  agg <- agg[order(agg$importance,decreasing = TRUE),]
  agg <- agg[grepl ("eig", agg$var),]
  return(agg$var[which.min(agg$importance)])}

#least_imp but using default hyper parameters and only considering ecomorph vars
least_imp_default_eco <- function(dat, wgt) {
  rf <- ranger(formula = interact ~ ., 
               data = dat,
               num.threads = 20,
               case.weights = wgt,
               probability = T,
               splitrule = "gini", #not much difference, but "gini' appears to do a bit better than "hellinger"
               importance = 'permutation')
  imp <- data.frame(importance(rf)) 
  #remove "source" and "target" so can aggregate importance
  imp$var<- gsub("source","", gsub("target", "", rownames(imp)))
  #aggregate
  agg <- aggregate(imp$importance.rf.,by=list(imp$var),FUN=mean)
  names(agg) <- c("var","importance")
  agg <- agg[order(agg$importance,decreasing = TRUE),]
  agg <- agg[!(grepl ("eig", agg$var)),]
  return(agg$var[which.min(agg$importance)])}
