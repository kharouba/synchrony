rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")
library(ggplot2)
library(rstan)
library(shinyStan)
best <- read.csv("stan_synchro.csv", header=TRUE)

studs <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/studies.csv", header=TRUE)
best2<-merge(best, unique(studs[,c("studyid","terrestrial2")]), by=c("studyid"))
#best2$terrestrial2<-as.factor(best2$terrestrial2)
best2$diff<-as.numeric(best2$diff)
best<-best2

N <- nrow(best)
K<-1
y <- best$diff
length <- best$length
la<- model.matrix(~ terrestrial, data = best)
eco<-la[,2]; eco<-as.vector(eco)
best[best$terrestrial2=="terrestrial",11]<-"0"
best[best$terrestrial2=="marine",11]<-"1"
best[best$terrestrial2=="fresh",11]<-"2"
best$terrestrial2<-as.integer(best$terrestrial2)
eco2<-as.vector(best$terrestrial2)
#X<-model.matrix(diff~terrestrial, data=best)

nVars<-1
Imat <- diag(1, nVars)
fits <- stan("quant_synch.stan", data=c("N","y","length","nVars","Imat"), iter=2000, chains=4)
fits <- stan("quant_synch.stan", data=c("N","y","eco2","nVars","Imat"), iter=2000, chains=4)
#fits <- stan("quant_synch.stan", data=c("N","K","y","X"), iter=2000, chains=4)
zz<-summary(fits)
launch_shinystan(fits)


