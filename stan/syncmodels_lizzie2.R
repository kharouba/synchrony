### Started 1 April 2015 ###
### April Fool's Day! ###

## Trying to run STAN with synchrony data ##

options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")
library(ggplot2)
library(rstan)
library(shinyStan)
# library(lme4)

# in the data file: spp1 = neg= negative species e.g. resource 
# data <- read.csv("input/alldata_lizzie.csv", header=TRUE)
# studies <- read.csv("input/studies.csv", header=TRUE)

# clean raw data
# in the data file: spp1 = neg= negative species e.g. resource 
raw <- read.csv("raw.csv", header=TRUE)
rawlong <- read.csv("indiv_sppdata2.csv", header=TRUE)
rawlong$int_type[which(rawlong$int_type=="poll")] <- "pollination"
rawlong.X <- NULL


##
# figure out how many unique species
# and alter the data so each unique species shows up once
# watch out: if you leave in int_type you still get duplicated species...
# so it's not below until after duplicate species are removed!
specieschar.wdups <- aggregate(rawlong["phenovalue"],
    rawlong[c("studyid", "species", "spp")], FUN=length)
specieschar <- aggregate(specieschar.wdups["phenovalue"],
    specieschar.wdups[c("studyid", "species")], FUN=length)
dupspp <- subset(specieschar, phenovalue>1)
specieschar.wdups[which(specieschar.wdups$species %in% dupspp$species),] #there are some species with mult relationships
# delete duplicate species as spp1 (generally the shorter timeseries)
rawlong.nodups <- rawlong[-(which(rawlong$species %in% dupspp$species &
    rawlong$spp=="spp1")),]
# and order it!
rawlong.nodups <- rawlong.nodups[with(rawlong.nodups, order(species, year)),]
specieschar.formodel <- aggregate(rawlong.nodups["phenovalue"],
    rawlong.nodups[c("studyid", "species", "int_type", "spp")], FUN=length)

# lme4: fails to converge
#modresource <- lmer(neg_phenovalue ~ year + (1 + year|spp1), data=raw)
#modresource <- lmer(neg_phenovalue ~ year + (year|spp1), data=raw)

##
# get estimate from no pooling (fixed effects) and complete pooling
comppool <- lm(phenovalue~year, data=rawlong.nodups)
comppool.wtype <- lm(phenovalue~year*spp-1, data=rawlong.nodups)
nopool <- lm(phenovalue~year*species, data=rawlong.nodups)
uniquespp <- unique(rawlong.nodups$species)
lmfits <- rep(NA, length(uniquespp))
for (eachsp in 1:length(uniquespp)){
    lmhere <- lm(phenovalue~year, data=subset(rawlong.nodups, species==uniquespp[eachsp]))
    lmfits[eachsp] <- coef(lmhere)[2]
  }



# try Stan
# 3 stages of cooking: prepared foods, recipes, make your own
# Stan does not currently allow ragged arrays
# see files Andrew sent on 1 April 2015
# mu is level as in level vs trend

# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues
rawlong.nodups$yr1976 <- rawlong.nodups$year-1976

N <- nrow(rawlong.nodups)
y <- rawlong.nodups$phenovalue
J <- nrow(specieschar.formodel)
species <- as.numeric(as.factor(rawlong.nodups$species))
year <- rawlong.nodups$yr1976
#type <- as.numeric(as.factor(specieschar.formodel$spp))

#Run the model
fit1<-stan("synchrony1_2.stan",data=c("N","y","J","species","year"), iter=100, chains=4)
print(fit1) #Rhat should be ~1

#level=intercept
#trend=slope
#sigma=variance

launch_shinystan(fit1)

#Hack to look at the data
qq<-summary(fit1) #gives the overall summary and then each of the four chains
qq[[1]][,1] #gives overall summary
spptrends.fromstan<-as.vector(qq[[1]][,1])[2:72]
tr<-cbind(lmfits, spptrends.fromstan); tr<-as.data.frame(tr)
spptrends$doybydec<-spptrends$spptrends.fromstan*10

ggplot(data=spptrends, aes(x=variable, y=doybydec))+
geom_path(aes(group=factor(intid)))+
#geom_point(aes(colour=factor(interaction), size=1))
theme_bw()+theme(panel.grid.major=element_blank(), axis.title.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))

## extra ##
# Plot the data
pdf("realdata.pdf", height=4, width=6)
colors=c("blue", "red")
par(mar=c(3,3,1,1), mgp=c(1.5,.5,0), tck=-.01)
plot(neg_phenovalue~year, data=raw, type="n", xlab="Year", ylab="Day of year", bty="l", main="Raw real data (resources and consumers are blue and red)")
for (j in 1:length(unique(raw$spp1))){
  subbydata <- subset(raw, spp1==unique(raw$spp1)[j])
  lines(neg_phenovalue~year, data=subbydata, col=colors[1])
}
for (j in 1:length(unique(raw$spp2))){
  subbydata <- subset(raw, spp2==unique(raw$spp2)[j])
  lines(pos_phenovalue~year, data=subbydata, col=colors[2])
}
dev.off()