### Plotting synchange versus tempchange ###
### Quick plot started 23 September 2016 ###

rm(list=ls())
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)

# setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016")

syncdat.full <- read.csv("output/sync.change.1K.csv", header=TRUE)
tempdat.full <- read.csv("output/temp.change.1K.csv", header=TRUE, col.names=c("rowhere","intid", "iteration", "temp.change"))

syncdat <- subset(syncdat.full, select=c("intid", "sync.change", "iteration"))
tempdat <- subset(tempdat.full, select=c("intid", "temp.change", "iteration"))

tempdat$iteration <- as.numeric(as.factor(tempdat$iteration))
mode(syncdat$iteration)

tryme <- inner_join(syncdat, tempdat, by=c("intid", "iteration"))
