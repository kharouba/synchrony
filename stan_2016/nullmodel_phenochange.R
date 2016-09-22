#Aug 2016- Look at workflow from Lizzie- 2016Aug22_nullmodel.pdf
## Null model used for individual species' phenological change


#Step 1- create pre_climate change dataset

#load rawlong.tot
rawlong.tot$count<-1
pre<-subset(rawlong.tot, year<=1981)
sss<- aggregate(pre["count"], pre[c("studyid", "intid", "species")], FUN=sum)
sss2<-subset(sss, count>=5) #datasets with enough data pre climate change
sss2$speciesid<-1:nrow(sss2) #number datas

pre_cc<-merge(rawlong.tot, sss2, by=c("studyid", "intid", "species"))

#Step 2- Create distribution of means and sd (2)
means <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "intid", "species")], FUN=mean)
names(means)[4]<-"mean_doy" #mean of each dataset
sds <- aggregate(pre_cc["phenovalue"], pre_cc[c("studyid", "intid", "species")], FUN=sd)
names(sds)[4]<-"sd_doy" #sd of each dataset

# FOR EACH species
clean<-unique(rawlong.tot[,c("studyid","int_type","year","species","phenofreq","short_site","terrestrial","phenovalue","newyear","yr1981","count")])
Bgroups<-unique(clean$species); b<-Bgroups; b<-as.character(b)
new<-data.frame(array(0,c(nrow(clean),12)))
rowcount<-1
for(i in 1:length(b)){
spp_long<-subset(clean, species==b[i])
asdf<-rowcount+(nrow(spp_long)-1)

#Step 3- Create new distribution for species from sample mean and sd
mu_spp1<-sample(means$mean_doy, size=1) #draw one mean from means
sigma_spp1<-sample(sds$sd_doy, size=1) #draw one sd from distribution of SDs
a<-rnorm(1000, mu_spp1, sigma_spp1); #create distribution based on moean and sd

#Step 3- Create new time-series for species of X length based on normal distribution
spp_long$ynull<-sample(a, size=nrow(spp_long)); 

new[rowcount:asdf,]<-spp_long
rowcount<-rowcount+nrow(spp_long)
asdf<-rowcount+(nrow(spp_long)-1)
}
names(new)[12]<-"ynull"
names(new)[1:11]<-names(clean)[1:11]

# Step 5- Run Stan
N <- nrow(new)
y <- new$ynull
specieschar.hin<- aggregate(new["ynull"], new[c("studyid", "species", "int_type", "terrestrial")], FUN=length) #number of years per species
specieschar.hin <- specieschar.hin[with(specieschar.hin, order(species)),]
Nspp <- nrow(specieschar.hin)
species <- as.numeric(as.factor(new$species))
year <- new$yr1981

null.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)

goo <- extract(null.model) 
specieschar.formodel.sm <- subset(specieschar.formodel, select=c("studyid", "species"))
intid <- read.csv("input/raw_april.csv", header=TRUE)
lal<-unique(rawlong.tot[,c("intid","terrestrial")])
intid2<-merge(intid, lal, by=c("intid"))
intid.sm <- subset(intid2, select=c("studyid", "spp1", "spp2", "intid" , "interaction","terrestrial"))
intid.nodups <- intid.sm[!duplicated(intid.sm),]
sync_int<-intid.nodups #synchrony change interactions
sync_int <- sync_int[with(sync_int, order(species)),]


summ_studyspp <- subset(specieschar.formodel, select=c("studyid", "species")); summ_studyspp<-unique(summ_studyspp)

## To estimate pheno change per species
it1000 <- matrix(0, ncol=3000, nrow=Nspp)
for (i in 3000:6000){ # 3000 iterations?
    summ_studyspp$model <- goo$b[i,]
    it1000[,(i-3000)] <- goo$b[i,]
}
summ_studyspp$stanfit <- rowMeans(it1000, na.rm=TRUE) #mean across iterations for EACH SPP


mean(summ_studyspp$stanfit) 

#computation of the standard error of the mean
sem<-sd(summ_studyspp$stanfit)/sqrt(length(summ_studyspp$stanfit)); sem
#95% confidence intervals of the mean
c(mean(summ_studyspp$stanfit)-2*sem,mean(summ_studyspp$stanfit)+2*sem)

