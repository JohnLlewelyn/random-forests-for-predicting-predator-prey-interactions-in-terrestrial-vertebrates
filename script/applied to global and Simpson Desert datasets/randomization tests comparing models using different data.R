#randomization tests comparing random forest performance predicting predator-prey interactions
#set file paths to run (search for lines with "###")

library(simstudy)
library(data.table)

#get data on performance of optimised models applied to global and Simpson Desert data 
# Make a vector of all your file paths
file_paths <- list.files(path = here::here("~/###/#WHERE RESULTS FROM APPLIED OPTIMISED MODELS SAVED#/"), pattern = "\\.rds", full.names = TRUE)
# Make a vector of file names
file_names <-  gsub(pattern = "\\.rds$", replacement = "", x = basename(file_paths))
# Read all your data into a list
data_list <- lapply(file_paths, readRDS)
# Assign file names to list elements
names(data_list) <- file_names   

#save as separate objects
list2env(data_list,envir=.GlobalEnv)

#reformat Simpson Desert data
SDaspplied_Eco10 <- unlist(SDaspplied_Eco10[1,])
SDaspplied_EcoAll <- unlist(SDaspplied_EcoAll[1,])
SDaspplied_PEMs10 <- unlist(SDaspplied_PEMs10[1,])
SDaspplied_PEMs10_Eco10 <- unlist(SDaspplied_PEMs10_Eco10[1,])
SDaspplied_PEMs21 <- unlist(SDaspplied_PEMs21[1,])
SDaspplied_PEMs21_EcoAll <- unlist(SDaspplied_PEMs21_EcoAll[1,])

#adapted from https://alfurka.github.io/2020-12-16-sampling-with-density-function/
#N.draw = 100000
samp_func <- function (x, N.draw){
dens = density(x)
cdf = cumsum(dens$y) / cumsum(dens$y)[length(dens$y)] #a random draw from the cumulative distribution function (CDF) that finds the index of the draw. I add 1 to the index because it finds the index value on the left. From ##dens$x[] I find the value of the drawn index.
samp = replicate(N.draw, dens$x[findInterval(runif(1), cdf)+1])
return(samp)}

#Comparing models with different variables
#GloBI trained and GloBI tested
gE10samp <- samp_func(GloBIonly_EcoA10,100000)
gP10samp <- samp_func(GloBIonly_PEMs10,100000)
gB20samp <- samp_func(GloBIonly_PEMs10_Eco10,100000)

dt <- data.frame(cbind(gE10samp,gP10samp,gB20samp))
dt$ecoVSpem <- dt$gE10samp>=dt$gP10samp
dt$ecoVSboth <- dt$gE10samp>=dt$gB20samp

table(dt$ecoVSpem)
table(dt$ecoVSboth)

#Comparing models with different variables
#GloBI trained and Simpson Desert tested
sE10samp <- samp_func(SDaspplied_Eco10,100000)
sP10samp <- samp_func(SDaspplied_PEMs10,100000)
sB20samp <- samp_func(SDaspplied_PEMs10_Eco10,100000)

dt <- data.frame(cbind(sE10samp,sP10samp,sB20samp))
dt$ecoVSpem <- dt$sE10samp>=dt$sP10samp
dt$ecoVSboth <- dt$sE10samp>=dt$sB20samp

table(dt$ecoVSpem)
table(dt$ecoVSboth)

sE21samp <- samp_func(SDaspplied_EcoAll,100000)
sP21samp <- samp_func(SDaspplied_PEMs21,100000)
sB42samp <- samp_func(SDaspplied_PEMs21_EcoAll,100000)

dt <- data.frame(cbind(sE21samp,sP21samp,sB42samp))
dt$ecoVSpem <- dt$sE21samp>=dt$sP21samp
dt$ecoVSboth <- dt$sE21samp>=dt$sB42samp

table(dt$ecoVSpem)
table(dt$ecoVSboth)

#compare many versus few ecomorphological models
gE10samp <- samp_func(GloBIonly_EcoA10,100000)
gE21samp <- samp_func(GloBIonly_EcoAll,100000)
table(gE10samp>=gE21samp)

sE10samp <- samp_func(SDaspplied_Eco10,100000)
sE21samp <- samp_func(SDaspplied_EcoAll,100000)
table(sE10samp>=sE21samp)

#################################################################
#compine the GloBI "few variable" model outputs into two separate df with groups to compare
Eco10 <- data.frame(tss=GloBIonly_EcoA10,grp = "Eco10")
PEM10 <- data.frame(tss=GloBIonly_PEMs10, grp = "PEM10")
both10 <- data.frame(tss=GloBIonly_PEMs10_Eco10, grp = "bothFew")
ecoVSpem <- rbind(Eco10,PEM10)
ecoVSboth <- rbind(Eco10, both10)
pVSb <- rbind(PEM10,both10)
colnames(ecoVSpem)<-c("Y","rx")
colnames(ecoVSboth)<-c("Y","rx")
colnames(pVSb)<-c("Y","rx")



