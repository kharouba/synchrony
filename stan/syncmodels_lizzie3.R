### Started 1 April 2015 ###
### April Fool's Day! ###

## Trying to run STAN with synchrony data ##

rm(list=ls()) 
options(stringsAsFactors = FALSE)

setwd("/users/kharouba/google drive/UBC/synchrony project/analysis/stan")
library(ggplot2)
library(rstan)
library(shinyStan)
library(reshape)
# library(lme4)

# in the data file: spp1 = neg= negative species e.g. resource 
# data <- read.csv("input/alldata_lizzie.csv", header=TRUE)
# studies <- read.csv("input/studies.csv", header=TRUE)

### clean raw data #####
# in the data file: spp1 = neg= negative species e.g. resource 
raw <- read.csv("raw_aug.csv", header=TRUE)
rawlong <- read.csv("indiv_sppdata_aug.csv", header=TRUE)
rawlong <- rawlong[with(rawlong, order(intid)),]
rawlong$int_type[which(rawlong$int_type=="poll")] <- "pollination"
rawlong.X <- NULL
change<-subset(rawlong, species=="Parus2 caeruleus" & intid=="216")
change$species[which(change$species=="Parus2 caeruleus")] <- "Parus2a caeruleus"

change2<-subset(rawlong, species=="Parus2 caeruleus" & intid=="220")
change2$species[which(change2$species=="Parus2 caeruleus")] <- "Parus2b caeruleus"

ori<-subset(rawlong, species!="Parus2 caeruleus")
ori2<-rbind(ori, change, change2)
rawlong<-ori2

change<-subset(rawlong, species=="Parus4 major" & intid=="221")
change$species[which(change$species=="Parus4 major")] <- "Parus4a major"

change2<-subset(rawlong, species=="Parus4 major" & intid=="217")
change2$species[which(change2$species=="Parus4 major")] <- "Parus4b major"

ori<-subset(rawlong, species!="Parus4 major")
ori2<-rbind(ori, change, change2)
rawlong<-ori2

#notes
#caterpillar spp.
#int177: 1961-2004
#int178: 1962-2006
#int179: 1961-2007
#for int209- Pleurobrachia pileus: 1976-2004 vs int208- 1976-2003
#Ficedula hypoleuca

# changed spp name because slightly different time series than other interaction in ms, now matches its interacting spp.
rawlong$rowid<-1:nrow(rawlong)
change<-subset(rawlong, species=="Caterpillar spp." & intid=="177")
change$species[which(change$species=="Caterpillar spp.")] <- "Caterpillar_a spp."
change$rowid<-2236:(2236+nrow(change)-1)
ori<-subset(rawlong, rowid<479 | rowid>514)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar spp." & intid=="178")
change$species[which(change$species=="Caterpillar spp.")] <- "Caterpillar_b spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<551 | rowid>584)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar spp." & intid=="179")
change$species[which(change$species=="Caterpillar spp.")] <- "Caterpillar_c spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<653 | rowid>686)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Pleurobrachia pileus" & intid=="209")
change$species[which(change$species=="Pleurobrachia pileus")] <- "Pleurobrachia_a pileus"
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1623 | rowid>1649)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Ficedula hypoleuca" & intid=="214")
change$species[which(change$species=="Ficedula hypoleuca")] <- "Ficedula_a hypoleuca"
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1677 | rowid>1692)
ori2<-rbind(ori, change)
rawlong<-ori2

#not different length of time series but easier to fix interactions HERE
#for int215- Parus ater
change<-subset(rawlong, species=="Parus ater" & intid=="215")
change$species[which(change$species=="Parus ater")] <- "Parus_a ater"
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1709 | rowid>1724)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="218")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2a spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1773 | rowid>1791)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="219")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2b spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1811 | rowid>1830)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="220")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2c spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1851 | rowid>1870)
ori2<-rbind(ori, change)
rawlong<-ori2

change<-subset(rawlong, species=="Caterpillar2 spp." & intid=="221")
change$species[which(change$species=="Caterpillar2 spp.")] <- "Caterpillar2d spp."
change$rowid<-max(rawlong$rowid):(max(rawlong$rowid)+nrow(change)-1)
ori<-subset(rawlong, rowid<1871 | rowid>1890)
ori2<-rbind(ori, change)
rawlong<-ori2

#check interactions to see if there are still 2 spp per interaction
#rename Parus4 major of int221

sub<-subset(rawlong, studyid=="HMK018")
sub$phenovalue2<-with(sub, phenovalue+182) #number of days different between Dec 22 (SH summer solstice) and June 21 (NH summer solstice), to put all dates on same 'calendar'
sub$phenovalue3<-with(sub, phenovalue2-365) #now same number of days before NH summer solstice as it was before SH summer solstice

rawlong2<-subset(rawlong, studyid!="HMK018")
rawlong2$phenovalue2<-rawlong2$phenovalue
rawlong2$phenovalue3<-rawlong2$phenovalue
rawlong3<-rbind(rawlong2, sub)
rawlong<-rawlong3
rawlong<-rawlong[,c("studyid","intid","int_type","spp","year","species","phenofreq","short_site","terrestrial","phenovalue3")]
names(rawlong)[10]<-"phenovalue"


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


##### MODELS #####
## Question 1 ##
rawlong.nodups$yr1981 <- rawlong.nodups$year-1981

N <- nrow(rawlong.nodups)
y <- rawlong.nodups$phenovalue
J <- nrow(specieschar.formodel)
species <- as.numeric(as.factor(rawlong.nodups$species))
year <- rawlong.nodups$yr1981
#y_prior_mean<-183 #median between 1 and 365
#y_prior_sd<-0.5 # concentrates 95% of data in interval, pg 43 in stan reference manual
#type <- as.numeric(as.factor(specieschar.formodel$spp))

#Run the model

# Margaret Kosmala's model with no type
fit.notype <- stan("synchrony1_notype.stan", data=c("N","y","J","species","year"), iter=2000, chains=3)
print(fit.notype)

#with covariance structure
nVars<-1
Imat <- diag(1, nVars)
fit.notype <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year", "nVars","Imat"), iter=2000, chains=4)

## TO EXTRACT SLOPE VALUES FROM MODEL##
#BIGGEST LOOP- FOCUSES ON SINGLE ITERATION FOR ALL SPECIES
asdf<-as.data.frame(unique(rawlong.nodups$species))
asdf$sppid<-1:nrow(asdf); names(asdf)[1]<-"species"
asdf2<-merge(asdf, unique(rawlong.nodups[,c("species","studyid","intid","spp")]), by="species")
tt<-as.data.frame(fit.notype)
rowcount<-1
lol<-rowcount+(54-1)
all<-data.frame(array(0, c(4000*87,4)))
names(all)[1]<-"intid"; names(all)[2]<-"iteration"; names(all)[3]<-"studyid"; names(all)[4]<-"diff"
for(i in 1:5){
tt2<-tt[i,88:174] # tt[1,88:175]
tt3<-t(tt2)

# to get data in an easier format to work with
new<-data.frame(array(0, c(87,3)))
names(new)[1]<-"iteration"; names(new)[2]<-"sppid";  names(new)[3]<-"slope"
new[,1]<-i
new[,2]<-asdf$sppid
new[,3]<-tt3

# for each iteration, take slope differences between interacting species
rawlong.nodups$sppid<-species
fits2<-merge(new, rawlong.nodups, by="sppid")
fits3<-unique(fits2[,c("sppid","iteration", "slope","intid", "spp","species","studyid")])
fits3 <- fits3[with(fits3, order(intid, spp)),]
names(fits3)[3]<-"slope"
#some species in multiple interactions WITHIN study and currently missing- need to add to fits3
#for int169- Pygoscelis antarcticus
newrow169<-c(61,i, tt2[[61]], 169, "spp2", "Pygoscelis antarcticus", "HMK018")
#newrow169<-as.data.frame(newrow169); as<-t(newrow169); as169<-as.data.frame(as); names(as169)<-names(fits3)
newrow171<-c(23, i, tt2[[23]], 171, "spp1", "Daphnia3 spp.", "HMK019")
#newrow171<-as.data.frame(newrow171); as<-t(newrow171); as171<-as.data.frame(as);names(as171)<-names(fits3)  
fits4<-rbind(fits3, as169, as171)
#fits4 <- fits4[with(fits4, order(intid, spp)),]
#fits3<-fits4
fits3$slope<-as.numeric(fits3$slope)
best<-data.frame(array(0, c(length(unique(fits3$intid)), 4)))
#names(best)[1]<-"intid"; names(best)[2]<-"iteration"; names(best)[3]<-"studyid"; names(best)[4]<-"diff"
Bgroups<-unique(fits3$intid); b<-Bgroups; b<-as.character(b)
for(j in 1:length(b)){
spp<-fits3[fits3$intid==b[j],]
best[j,1]<-b[j]
best[j,2]<-spp[1,2]
best[j,3]<-spp[1,7]
best[j,4]<-spp[1,c("slope")]-spp[2,c("slope")] #resource-consumer
}
#cut and paste each iteration and slope differences into large datasheet with all other iterations
all[rowcount:lol,]<-best
rowcount<-rowcount+(nrow(best)-1)
lol<-rowcount+(nrow(best)-1)
}

log_lik<-extract_log_lik(fit.notype, parameter_name = "log_lik")

colVars <- function (a){
  diff <- a - matrix (colMeans(a), nrow(a), ncol(a), byrow=TRUE)
  vars <- colMeans (diff^2)*nrow(a)/(nrow(a)-1)
  return (vars)
}

waic <- function (stanfit){ 
log_sum_exp <- function(x) { 
	x_max <- max(x) 
	x_max + log(sum(exp(x - x_max))) 
}  
log_lik <- extract (stanfit, "log_lik")$log_lik 
summands <- apply(log_lik, 2, log_sum_exp) 
correc <- - ncol(log_lik) * log(nrow(log_lik)) 
lppd <- sum(summands) + correc 
p_waic_1 <- 2 * (sum(summands - colMeans(log_lik)) + correc) 
p_waic_2 <- sum (colVars(log_lik)) 
waic_2 <- -2*lppd + 2*p_waic_2 
return (list (waic=waic_2, p_waic=p_waic_2, lppd=lppd, p_waic_1=p_waic_1)) 
} 
print (waic (fit.notype))


#alpha=intercept
#beta=slope
zz<-summary(fit.notype)
intercepts <- as.vector(zz[[1]][,1])[1:87]
spptrends<- as.vector(zz[[1]][,1])[88:174]
sppSE<- as.vector(zz[[1]][,2])[88:174]
sppid<-1:87
rawlong.nodups$sppid<-species
#fits<-data.frame(array(0, c(length(unique(rawlong.nodups$sppid)+1), 2)))
fits<-data.frame(array(0, c(length(unique(rawlong.nodups$sppid)+1), 2)))
fits$sppid<-sppid; fits$spptrends<-spptrends; fits$sppSE<-sppSE; fits<-fits[,3:5]
fits2<-merge(fits, rawlong.nodups, by="sppid")
fits3<-unique(fits2[,c("sppid","spptrends","sppSE","intid", "spp","species","studyid","int_type","phenofreq","short_site","terrestrial")])
fits3 <- fits3[with(fits3, order(intid, spp)),]

#some species in multiple interactions WITHIN study and currently missing- need to add to fits3
#for int169- Pygoscelis antarcticus
newrow169<-c(61,-0.2640256973, 0.009129462, 169, "spp2", "Pygoscelis antarcticus", "HMK018", "competition", "daily", "antarctica", "terrestrial")

#search for spp when included in another interaction (but same paper) to get slope but use info from current interaction to get other info
sub<-subset(fits3, species=="Daphnia3 spp.")
newrow171<-c(23, -0.6614358, 0.002896469, 171, "spp1", "Daphnia3 spp.", "HMK019", "predation","daily","windermere","aquatic")


fits4<-rbind(fits3, newrow169, newrow171)
fits4 <- fits4[with(fits4, order(intid, spp)),]
fits3<-fits4
fits3$spptrends<-as.numeric(fits3$spptrends)

# POS diff (diff>0)= decreasing synchrony; NEG diff (diff<0)= # increasing synchrony

best<-data.frame(array(0, c(length(unique(fits3$intid)), 7)))
names(best)[1]<-"intid"; names(best)[2]<-"studyid"; names(best)[3]<-"int_type"; names(best)[4]<-"phenofreq"; names(best)[5]<-"short_site"; names(best)[6]<-"terrestrial"; names(best)[7]<-"diff"
Bgroups<-unique(fits3$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-fits3[fits3$intid==b[i],]
best[i,1]<-b[i]
best[i,2:6]<-spp[1,7:11]
best[i,7]<-spp[1,c("spptrends")]-spp[2,c("spptrends")] #resource-consumer
}

obs<-mean(best$diff, na.rm=T)
new <- new[with(new, order(null_diff)),]
new$rank<-1:nrow(new)
obs has rank position of 982/983

# overall model #
sums<-merge(best, unique(raw[,c("intid","length","minyear")]), by=c("intid"))
newrow169<-c("221", "HMK048","predation", "weekly","netherlands", "terrestrial", -0.21475846, 20, 1985)
newrow170<-c("222", "HMK048","herbivory", "weekly","netherlands", "terrestrial", 0.07032723, 17, 1988)
sums2<-rbind(sums, newrow169, newrow170)
sums2$length<-as.integer(sums2$length)
sums2$minyear<-as.integer(sums2$minyear)
write.csv(sums2, "stan_synchro.csv")
m1<-lm(diff~length, data=sums2); summary(m1)
m1<-lm(diff~minyear, data=sums2); summary(m1)
m1<-lm(diff~terrestrial, data=sums2); summary(m1)
studs <- read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/studies.csv", header=TRUE)
sums3<-merge(sums2, unique(studs[,c("studyid","terrestrial2")]), by=c("studyid"))
sums3$terrestrial2<-as.factor(sums3$terrestrial2)
sums3$diff<-as.numeric(sums3$diff)
m1<-lm(diff~terrestrial2, data=sums3); summary(m1)

library(arm)
m1<-bayesglm(spptrends~terrestrial, family=gaussian, prior.mean=0, prior.scale=Inf, prior.df=Inf, data=fits3)
simulates<-coef(sim(m1))
post.open<-simulates[,2]
quantile(post.open, c(0.025,0.975)) # to get 95 credible intervals
library(MCMCpack) #need to get MCMC simulation method of creating emperical distribution
m2<-MCMCregress(spptrends~terrestrial, burnin=3000, mcmc=10000, thin=1, verbose=0, seed=NA, beta.start=NA, data=fits3)


#print(fit1) #Rhat should be ~1

#level=intercept
#trend=slope
#sigma=variance

# Predictive checks #
pred.fit <- stan("predcheck.stan", data=c("N","y","J","species","year"), iter=2000, chains=3)

n.dogs <- J;
n.trials <- 37 #(max number of years per species across interactions)
data <- list ("y", "n.dogs", "n.trials")
inits <- function (){list (b.0=rnorm(1), b.1=rnorm(1), b.2=rnorm(1)) } 
parameters <- c ("b.0", "b.1", "b.2", "y.rep", "y.mean.rep")
fit.logit.1 <- bugs(data, inits, parameters, "dogs.logit.1.bug", n.chains=3, n.iter=2000, n.thin=100) 

n.sims<-30
y.rep <- array (NA, c(n.sims, n.dogs, n.trials)) 
for (j in 1:n.dogs){
n.avoid.rep <- rep (0, n.sims)
n.shock.rep <- rep (0, n.sims) 
for (t in 1:n.trials){ 
p.rep <- invlogit (b.0 + b.1*n.avoid.rep + b.2*n.shock.rep)
y.rep[,j,t] <- rbinom (n.sims, 1, p.rep)
n.avoid.rep <- n.shock.rep + 1 - y.rep[,j,t]
n.shock.rep <- n.shock.rep + y.rep[,j,t] 
} } 



launch_shinystan(fit.nohinge)

## Question 2- is temp change driving changes in pheno synchrony? ##
rawlong.nodups2<-subset(rawlong.nodups, studyid!="HMK029" | studyid!="HMK037" | studyid!="HMK032")
rawlong.nodups<-rawlong.nodups2

#nohinge
rawlong.nodups$yr1981 <- rawlong.nodups$year-1981
N <- nrow(rawlong.nodups)
y <- rawlong.nodups$phenovalue
specieschar.no<- aggregate(rawlong.nodups["phenovalue"], rawlong.nodups[c("studyid", "species", "int_type", "spp")], FUN=length)
J <- nrow(specieschar.no)
species <- as.numeric(as.factor(rawlong.nodups$species))
year <- rawlong.nodups$yr1981
nVars<-1
Imat <- diag(1, nVars)
fit.nohinge <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=3)



#hinge model
#figure out how many interactions have n=>5 yrs before 1981
rawlong.nodups2<-subset(rawlong.nodups, year<=1981)
rawlong.nodups2$count<-1 # if need to count, subset to spp2
wus <- aggregate(rawlong.nodups2["count"], rawlong.nodups2[c("intid")], FUN=sum)

hinge<-subset(rawlong.nodups, intid=="170" | intid=="171" | intid=="177" | intid=="178" | intid=="179" | intid=="180" |intid=="181" | intid=="189" | intid=="191"|intid=="193" |intid=="194" | intid=="195" |intid=="196" |intid=="201" |intid=="207" |intid=="208")
hinge_non<-subset(rawlong.nodups, intid!="170" & intid!="171" & intid!="177" & intid!="178" & intid!="179" & intid!="180" & intid!="181" & intid!="189" & intid!="191" & intid!="193" & intid!="194" & intid!="195" & intid!="196" & intid!="201" & intid!="207" & intid!="208")
hinge_non$newyear<-hinge_non$year
hinge_pre<-subset(hinge, year<=1981); hinge_pre$newyear<-1981
hinge_post<-subset(hinge, year>1981); hinge_post$newyear<-hinge_post$year
hinges<-rbind(hinge_pre, hinge_post)
rawlong.tot<-rbind(hinge_non, hinges)

rawlong.tot$yr1981 <- rawlong.tot$newyear-1981
N <- nrow(rawlong.tot)
y <- rawlong.tot$phenovalue
specieschar.hin<- aggregate(rawlong.tot["phenovalue"], rawlong.tot[c("studyid", "species", "int_type", "spp")], FUN=length)
J <- nrow(specieschar.hin)
species <- as.numeric(as.factor(rawlong.tot$species))
year <- rawlong.tot$yr1981
nVars<-1
Imat <- diag(1, nVars)
fit.hinge <- stan("synchrony1_notype_hk.stan", data=c("N","J","y","species","year","nVars","Imat"), iter=2000, chains=4)

zzhinge<-extract(fit.hinge)
zzno<-extract(fit.nohinge)

zzhinge<-summary(fit.hinge)
intercepts <- as.vector(zzhinge[[1]][,1])[1:87]
spptrends<- as.vector(zzhinge[[1]][,1])[88:174]
sppSE<- as.vector(zzhinge[[1]][,2])[88:174]
sppid<-1:87
rawlong.nodups$sppid<-species
#fits<-data.frame(array(0, c(length(unique(rawlong.nodups$sppid)+1), 2)))
fits<-data.frame(array(0, c(length(unique(rawlong.nodups$sppid)+1), 2)))
fits$sppid<-sppid; fits$spptrends<-spptrends; fits$sppSE<-sppSE; fits<-fits[,3:5]
fits2<-merge(fits, rawlong.nodups, by="sppid")
fits3<-unique(fits2[,c("sppid","spptrends","sppSE","intid", "spp","species","studyid","int_type","phenofreq","short_site","terrestrial")])
fits3 <- fits3[with(fits3, order(intid, spp)),]



##log ratio:
(mean(zzhinge$lp__))/(mean(zzno$lp__))

posterior<- extract(fit.nohinge); 
mean(posterior$residual)
sigma.est <- mean(extract(fit.hinge, pars = 'sigma_b', inc_warmup = FALSE))
residual.est <- (extract(fit.hinge, pars = 'residual', inc_warmup = FALSE))
predicted.est <- (extract(fit.nohinge, pars = 'predicted', inc_warmup = FALSE))

#Hack to look at the data
qq<-summary(fit1) #gives the overall summary and then each of the four chains
qq[[1]][,1] #gives overall summary
#spptrends.fromstan<-as.vector(qq[[1]][,1])[2:72]
spptrends.fromstan<-as.vector(qq[[1]][,1])[72:142]
tr<-cbind(lmfits, spptrends.fromstan); tr<-as.data.frame(tr)
spptrends$doybydec<-spptrends$spptrends.fromstan*10

ggplot(data=spptrends, aes(x=variable, y=doybydec))+
geom_path(aes(group=factor(intid)))+
#geom_point(aes(colour=factor(interaction), size=1))
theme_bw()+theme(panel.grid.major=element_blank(), axis.title.x=element_text(size=12), axis.text.y=element_text(size=12), axis.title.y=element_text(size=12, angle=90))

## extra ##
spp1<-subset(fits3, spp=="spp1"); spp1<-unique(spp1[,c("sppid","spptrends","sppSE", "intid")]); names(spp1)[1]<-"sppid1"; names(spp1)[2]<-"spptrends1"; names(spp1)[3]<-"sppSE1"
spp2<-subset(fits3, spp=="spp2"); spp2<-unique(spp2[,c("sppid","spptrends","sppSE", "intid")])
tot<-merge(spp1, spp2, by=c("intid"))
tot$diff<-with(tot, spptrends1-spptrends)

dec:
dec1<-subset(tot, spptrends1>0 & spptrends>0 & diff>0)
dec2<-subset(tot, spptrends1<0 & spptrends<0 & diff>0)
dec3<-subset(tot, spptrends1>0 & spptrends<0)
dec<-rbind(dec1,dec2,dec3); nrow(dec)
inc:
inc1<-subset(tot, spptrends1>0 & spptrends>0 & diff<0)
inc2<-subset(tot, spptrends1<0 & spptrends<0 & diff<0)
inc3<-subset(tot, spptrends1<0 & spptrends>0)
inc<-rbind(inc1,inc2, inc3); nrow(inc)

##
# model with type (resource or consumer); the one and only original
fit1 <- stan("stan/synchrony1.stan", data=c("N","y","J","species","year","type"), iter=1000, chains=4)

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