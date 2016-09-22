# Sync models
rawlong.tot$yr1981 <- rawlong.tot$newyear-1981
uniquespp<-unique(rawlong.tot$species)
lmfits<-rep(NA< length(uniquespp))
ACfits<-rep(NA< length(uniquespp))
for(eachsp in 1:length(uniquespp)){
	lmhere<-gls(phenovalue~yr1981, data=subset(rawlong.tot, species==uniquespp[eachsp]))
	tempAC<-gls(phenovalue~yr1981, correlation=corAR1(), data=subset(rawlong.tot, species==uniquespp[eachsp])) 
	lmfits[eachsp]<-AIC(lmhere)
	ACfits[eachsp]<-AIC(tempAC)
	}
dr<-data.frame(cbind(lmfits, ACfits))
dr$lmfits<-as.numeric(dr$lmfits)
dr$ACfits<-as.numeric(dr$ACfits)
dr$AICdiff<-with(dr, lmfits-ACfits)

#Temp models
uniquespp<-unique(dataset$datasetid)
lmfits<-rep(NA< length(uniquespp))
ACfits<-rep(NA< length(uniquespp))
# No pooling:
for(eachsp in 1:length(uniquespp)){
	lmhere<-gls(envvalue~yr1981, data=subset(dataset, datasetid==uniquespp[eachsp]))
	tempAC<-gls(envvalue~yr1981, correlation=corAR1(), data=subset(dataset, datasetid==uniquespp[eachsp])) 
	lmfits[eachsp]<-AIC(lmhere)
	ACfits[eachsp]<-AIC(tempAC)
}
dr<-data.frame(cbind(lmfits, ACfits))
dr$lmfits<-as.numeric(dr$lmfits)
dr$ACfits<-as.numeric(dr$ACfits)
dr$AICdiff<-with(dr, lmfits-ACfits)

#Correlogram: ACF (only works with gls)
plot(ACF(tempAC), alpha=0.01) #looking for >1 lag to exceed the 95% significance bounds excluding first lag

sub<-subset(dataset, datasetid=="HMK011 _ 1")
#sub<-subset(dataset, datasetid=="HMK019 _ 3")
x<-sub$yr1981
y<-sub$envvalue
N<-nrow(sub)

lm.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/simplelinearmodel.stan", data=c("N","x","y"), iter=3000, chains=4)

alpha=-1.39
beta=0.15

lm<-lm(y~x); summary(lm)
sub$ypred<-predict(lm)

ar.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/ARmodel.stan", data=c("N","x","y"), iter=3000, chains=4)
asdf<-summary(ar.model, pars="beta")
median<-asdf[[1]][1:91];


ggplot(sub, aes(x=x, y=y))+geom_line()+geom_line(aes(x=x, y=ypred, colour="red"))+geom_abline(slope=0.32)

year<-dataset$yr1981
y<-dataset$envvalue
N<-nrow(dataset)
specieschar.hin<- aggregate(dataset["envvalue"], dataset[c("datasetid")], FUN=length) #number of years per dataset
specieschar.hin <- specieschar.hin[with(specieschar.hin, order(datasetid)),]
Nspp <- nrow(specieschar.hin) #J
species <- as.numeric(as.factor(dataset$datasetid))

AR.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2_AR.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)
temp.model<-stan("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/stanmodels/twolevelrandomslope2_AR.stan", data=c("N","Nspp","y","species","year"), iter=3000, chains=4)
