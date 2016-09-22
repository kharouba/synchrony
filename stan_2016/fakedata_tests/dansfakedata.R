nsite=2
nsp=20
nind=10
nwarm=2
nphoto=2
rep=1
ntot=nsite*nwarm*nphoto

#Build up data frame
site=gl(nsite, rep, length=ntot)
warm=gl(nwarm, rep*nsite, length=ntot)
photo=gl(nphoto, rep*nsite*nwarm, length=ntot)
treatcombo=paste(warm,photo, sep="_")
d<-data_frame(site, warm, photo, treatcombo)

# set up differences for each level
sitediff=2
warmdiff=-20
photodiff=-14
sitewarm=0
sitephoto=0
warmphoto=3.5

# interactions. 9 two-way interactions
sitediff.sd=1.5
warmdiff.sd=1
photodiff.sd=1
sitewarm.sd=1
sitephoto.sd=1
warmphoto.sd=1

mm<-model.matrix(~(site+warm+photo)^2, data.frame(warm, photo))

#with individuals
baseinter=35 #baseline intercept
spint<-baseinter+c(1:nsp)-mean(1:nsp) # diff intercepts by species

fake<-vector()

for(i in 1:nsp){ #loop over species, as these are the random effect modeled
	#within speices, have a loop for individuals
	indint<-spint[i]+1:nind-mean(1:nind) #diff intercept by individual
	
	for (j in 1:nind){
		coeff<-c(indint[j],
		rnorm(1,sitediff, sitediff.sd),
		rnorm(1, warmdiff, warmdiff.sd),
		rnorm(1, photodiff, photodiff.sd),
		rnorm(1, sitewarm, sitewarm.sd),
		rnorm(1, sitephoto, sitephoto.sd),
		rnorm(1, warmphoto, warmphoto.sd))
	bb<-rnorm(n=length(warm), mean=mm%*%coeff,sd=0.1)
	fakex<-data.frame(bb, sp=i, ind=paste(i,j,sep="_"), site, warm, photo)
	fake<-rbind(fake, fakex)
	}
}
