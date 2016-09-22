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
raw<-read.csv("raw.csv", header=TRUE)
Uria aalge}

rawlong<-read.csv("alldata_long.csv", header=TRUE)
rawlong<-rawlong[order(rawlong$value),]
new<-read.csv("indiv_sppdata.csv", header=TRUE)
Ficedula albicollis from 1,2
Mnemiopsis leidyi

# lme4: fails to converge
#modresource <- lmer(neg_phenovalue ~ year + (1 + year|spp1), data=raw)
#modresource <- lmer(neg_phenovalue ~ year + (year|spp1), data=raw)

# try Stan
# 3 stages of cooking: prepared foods, recipes, make your own
# Stan does not currently allow ragged arrays
# see files Andrew sent on 1 April 2015
# mu is level as in level vs trend

# prep the data to fit the model including:
# aggregate to get species level characteristics
# subset down to the phenovalues
longpheno<- subset(rawlong, whatmeasured=="phenovalue")
longpheno$yr1976 <-longpheno$year-1976
specieschar<-aggregate(longpheno["measurement"], longpheno[c("studyid","value","variable","intid")], FUN=length)
specieschar.nodup<-specieschar[!(duplicated(specieschar$value)==TRUE),] #there are some species with mult relationships

N<-nrow(longpheno)
y<-longpheno$measurement
J<-nrow(specieschar.nodup)
species<-as.numeric(as.factor(longpheno$value)) #converts spp names to numbers
year<-longpheno$yr1976
type<-as.numeric(as.factor(specieschar.nodup$variable)) # converts spp type (neg or pos) to 1 or 2

#Run the model
fit1<-stan("synchrony1.stan",data=c("N","y","J","species","year","type"), iter=100, chains=4)
print(fit1) #Rhat should be ~1

#level=intercept
#trend=slope
#sigma=variance

launch_shinystan(fit1)

#Hack to look at the data
qq<-summary(fit1) #gives the overall summary and then each of the four chains
qq[[1]][,1] #gives overall summary
spptrends.fromstan<-as.vector(qq[[1]][,1])[4:73]
spptrends<-cbind(specieschar.nodup, spptrends.fromstan)
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