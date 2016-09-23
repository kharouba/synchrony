### Plotting synchange versus tempchange ###
### Quick plot started 23 September 2016 ###

rm(list=ls())
options(stringsAsFactors=FALSE)

library(plyr)
library(dplyr)
library(ggplot2)

# setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016")
setwd("~/Documents/git/projects/trophsynch/synchrony/stan_2016")

syncdat.full <- read.csv("output/sync.change.1K.csv", header=TRUE)
tempdat.full <- read.csv("output/temp.change.1K.csv", header=TRUE, col.names=c("rowhere","intid", "iteration", "temp.change"))

syncdat <- subset(syncdat.full, select=c("intid", "sync.change", "iteration"))
tempdat <- subset(tempdat.full, select=c("intid", "temp.change", "iteration"))

tempdat$iteration <- as.numeric(as.factor(tempdat$iteration))
mode(syncdat$iteration)

compdat <- inner_join(syncdat, tempdat, by=c("intid", "iteration"))

## summarize the data (not the right way to do this, but at least quick and applicable)
datersumm <-
      ddply(compdat, c("intid"), summarise,
      mean.tc = mean(temp.change, na.rm=TRUE),
      sd.tc = sd(temp.change, na.rm=TRUE),
      mean.sc = mean(sync.change, na.rm=TRUE),
      sd.sc = sd(sync.change, na.rm=TRUE))

## get error bars quick and dirty
datersumm$max.tc <- datersumm$mean.tc + datersumm$sd.tc
datersumm$min.tc <- datersumm$mean.tc - datersumm$sd.tc
datersumm$max.sc <- datersumm$mean.sc + datersumm$sd.sc
datersumm$min.sc <- datersumm$mean.sc - datersumm$sd.sc

ggplot(datersumm, aes(x=mean.tc, y=mean.sc, colour=as.factor(intid))) +
    geom_point(size=1.25) +
    geom_errorbar(aes(ymin=min.sc, ymax=max.sc), size=0.75) + 
    geom_errorbarh(aes(xmin = min.tc,xmax = max.tc), size=0.75) + 
    xlab("temp change (C)") +
    ylab("synchrony change") +
    theme_bw() +
    theme(legend.key = element_blank())

ggplot(datersumm, aes(x=mean.tc, y=mean.sc, colour=as.factor(intid))) +
    geom_point(size=1.25)

plot(mean.sc~mean.tc, datersumm)
abline(lm(mean.sc~mean.tc, datersumm))
summary(lm(mean.sc~mean.tc, datersumm))

plot(abs(mean.sc)~mean.tc, datersumm)
abline(lm(abs(mean.sc)~mean.tc, datersumm))
summary(lm(abs(mean.sc)~mean.tc, datersumm))
