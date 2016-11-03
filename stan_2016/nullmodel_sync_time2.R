#Aug 2016- Look at workflow from Lizzie- 2016Aug22_nullmodel.pdf
## Null model used for synchrony change, to explore effect of length of time series on ability to detect sync change

Overall workflow
Create spp1 data based on pre climate change or all years of data
Run stan for spp1
Determine spp2 slope based on randomly sampled sync change from either observed distribution OR simulated distribution


#Step 1- create pre_climate change dataset
#load rawlong.tot
interact <- read.csv("input/raw_oct.csv", header=TRUE)
interact$newphenodiff<- with(interact, neg_phenovalue-pos_phenovalue) #spp1-spp2

rawlong.tot$count<-1
pre<-subset(rawlong.tot, year<=1981)
sums<-unique(rawlong.tot[,c("studyid","intid","species")])

#some species have non-unique data e.g. Glis glis, C

pre_uni<-unique(pre[,c("studyid","year","species","phenovalue","yr1981","count")])

#pres<-merge(pre, unique(rawlong.tot[,c("studyid","intid","species")]), by=c("studyid","species"))
sss<- aggregate(pre_uni["count"], pre_uni[c("studyid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(pre_uni[,1:5], sss2, by=c("studyid","species"))
#pre_cc2<-merge(pre_cc, unique(rawlong.tot[,c("studyid","intid","species")]), by=c("studyid","species"))
#pre_cc2 <- pre_cc2[with(pre_cc2, order(species, year)),]
#pre_cc<-merge(pres[,c(1:4,6)],  sss2, by=c("studyid", "intid", "species"))

OR Use all data
rawlong.tot2<-unique(rawlong.tot[,c("studyid","species","phenovalue","yr1981")]) #CLEAN UP so only unique values across repeating species within studoes
pre_cc<-rawlong.tot2


!!! Manual up to here!!


#Step 2- Create distribution of means and sd (2)
means <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=mean) #(not intid)
names(means)[3]<-"mean_doy" #mean of each dataset
sds <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "species")], FUN=sd)
names(sds)[3]<-"sd_doy" #sd of each dataset

# For each length of time series

timeseries<-c(5,10,15,20,25,30,35,40)
rowcount2<-1
asdf2<-rowcount2+(54-1) #number of interactions
output<-data.frame(array(0,c(8*54,4)))

for(le in 1:length(timeseries)){ #length(timeseries)


# FOR EACH INTERACTION
Bgroups<-unique(rawlong.tot$intid); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0,c(timeseries[le]*54,5)))
rowcount<-1
for(i in 1:length(b)){
int_long<-subset(rawlong.tot, intid==b[i])
asdf<-rowcount+(timeseries[le]-1) #length of NEW time series

#Step 3- Create new distribution for spp1 from sample mean and sd
mu_spp1<-sample(means$mean_doy, size=1) #draw one mean from means
sigma_spp1<-sample(sds$sd_doy, size=1) #draw one sd from distribution of SDs
a<-rnorm(1000, mu_spp1, sigma_spp1); #create distribution of phenovalues based on moean and sd


#Step 3- Create new time-series for spp1 of X length based on normal distribution
spp1_data<-subset(int_long, spp=="spp1")
time<-spp1_data$yr1981
time2<-rnorm(1000, mean(time), sd(time)) # create distribution of years based on observed
spp1_datab<-data.frame(array(0,c(timeseries[le],3))) #make new dataframe based on new length
spp1_datab[,1]<-rep(unique(spp1_data$studyid), timeseries[le])
spp1_datab[,2]<-rep(unique(spp1_data$intid), timeseries[le])
spp1_datab[,3]<-rep(unique(spp1_data$species), timeseries[le])

spp1_datab$ynull<-sample(a, size=timeseries[le]); #start with a 5 year time series
spp1_datab$yr1981<-sample(time2, size=timeseries[le]); #sample X number of years based on length of new time series

new[rowcount:asdf,]<-spp1_datab
rowcount<-rowcount+timeseries[le]
asdf<-rowcount+(timeseries[le]-1)
}
names(new)[1]<-"studyid"; names(new)[2]<-"intid"; names(new)[3]<-"species"; names(new)[4]<-"ynull"; names(new)[5]<-"yr1981"

#spa<-rbind(spp1_datab, spp2_datab)
#names(spa)[1]<-"studyid"; names(spa)[2]<-"intid"; names(spa)[3]<-"species"; names(spa)[4]<-"ynull"; names(spa)[5]<-"yr1981"
#ggplot(spa, aes(x=yr1981, y=ynull))+geom_point(aes(colour=factor(species)))+geom_smooth(method="lm", se=FALSE, aes(colour=factor(species)))

#Step 4- clean to match syncmodels code, gets rid of non-unique data for repeating species within a study (i.e. those species across multiple intxns, e.g. Accipiter nisus) so that models for unique species are built, not unique species-intxns ; Number of unique species should match syncmodels i.e. n=85
new <- new[with(new, order(species)),]
Bgroups<-unique(new$species); b<-Bgroups; b<-as.character(b)
evennewer<-data.frame(array(0,c(timeseries[le]*54,5)))
rowcount<-1
asdf<-rowcount+(timeseries[le]-1)
for(i in 1:length(b)){
spp<-subset(new, species==b[i])

if(length(unique(spp$intid))==1){
	spp<-subset(new, species==b[i])
	evennewer[rowcount:asdf,]<-spp
	rowcount<-rowcount+(timeseries[le])
	asdf<-rowcount+(timeseries[le]-1)
}
if(length(unique(spp$intid))>1){
	cgroups<-unique(spp$intid); c<-cgroups; c<-as.character(c)
	spp<-subset(new, species==b[i] & intid==c[sample(length(unique(spp$intid)), size=1)])
	evennewer[rowcount:asdf,]<-spp
rowcount<-rowcount+(timeseries[le])
asdf<-rowcount+(timeseries[le]-1)
}
}
names(evennewer)<-names(new)
evennewer2<-subset(evennewer, studyid!=0)

# Step 5- Run Stan

N <- nrow(evennewer2)
y <- evennewer2$ynull
Nspp <- length(unique(evennewer2$species))
J <- length(unique(evennewer2$species))
species <- as.numeric(as.factor(evennewer2$species))
year <- evennewer2$yr1981

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

yup.sim <- extract(null.model) 
uni<-unique(evennewer2[,c("studyid","species")])
uni$spp<-"spp1"

# For spp' phenological change
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    uni$model <- yup.sim$b[i,]
    it1000[,(i-3000)] <- yup.sim$b[i,]
}
uni$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP

# Step 6 calculate slope for spp2
all<-merge(unique(rawlong.tot[,c("studyid","intid","species","spp")]), uni, by=c("studyid","species","spp"), all.x=TRUE)
all  <- all[with(all , order(intid, spp)),]

# load andtheanswer from syncmodels i.e. go with unmanipulated estimates of sync change
# create distribution of sync change based on observed mean and sd synchange from data
obs_syncchange<-rnorm(1000, mean(andtheanswer$meanchange), sd(andtheanswer$meanchange)); #cr

Bgroups<-unique(all$intid); b<-Bgroups; b<-as.character(b)
for(i in 1:length(b)){
int_long<-subset(all, intid==b[i])
spp2<-int_long[1,"stanfit"]-sample(obs_syncchange, size=1) #draw one sync change from means AND subtract spp2 from spp1
small<-c(int_long[1,"stanfit"],sample(obs_syncchange, size=1)) #draw one sync change from means
all[all$intid==b[i],"stanfit"]<-small
}
