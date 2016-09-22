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
