total3<-read.csv("int_phenodata.csv", header=TRUE, na.strings="NA", as.is=TRUE)

### to deal with pseudo-replication for INTERACTIONS
#remove Keratella cochlearis from either HMK029 or HMK037
#total<-subset(total3, spp2!="Keratella1 cochlearis") # from HMK029
total<-subset(total3, intid!="182") #see right above for reason
total<-subset(total, studyid!="HMK026")
total<-subset(total, studyid!="HMK049") #HMK049- remove because predicted relationship


total3<-total
total3$intxn<-with(total3, paste(spp1,spp2,interaction))
total3<-subset(total3, intxn!="Caterpillar2 spp. Parus4 major predation") #same interaction as HMK027, excluded because shorter time series

total3<-subset(total3, intxn!="Quercus3 robur Caterpillar2 spp. herbivory")

#number of years by INTERACTION (not study!)
total4<-na.omit(total3)
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=length)
names(yo)[names(yo)=="year"]<-"length"
yo2<-merge(total4,yo, by=c("studyid","intid"))

#first year of study WITH INTERACTION
yo<-aggregate(total4["year"], by=total4[c("studyid","intid")], FUN=min)
names(yo)[names(yo)=="year"]<-"minyear"
yo3<-merge(yo2,yo, by=c("studyid","intid"))

#re-calculate phenodiff from min year
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
best<-data.frame(array(0, c(nrow(yo3), 4)))
names(best)[1]<-"intid"; names(best)[2]<-"year"; names(best)[3]<-"phenodiff_base"; names(best)[4]<-"base"
Bgroups<-unique(yo3$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-yo3[yo3$intid==b[i],]
asdf<-rowcount+(nrow(spp)-1)
best[rowcount:asdf,1]<-b[i]
best[rowcount:asdf,2]<-spp[,c("year")]
best[rowcount:asdf,3]<- with(spp, phenodiff-spp[1,c("phenodiff")])
best[rowcount:asdf,4]<- spp[1,c("phenodiff")]
rowcount<-rowcount+nrow(spp)
}
yo4<-merge(yo3, best, by=c("intid","year"))
yo3<-yo4

#re-calculate phenodiff from 3 years
yo3<-yo3[order(yo3$intid,yo3$year),]; yo3<-na.omit(yo3)
best<-data.frame(array(0, c(nrow(yo3), 3)))
names(best)[1]<-"intid"; names(best)[2]<-"year"; names(best)[3]<-"phenodiff_baseavg"
Bgroups<-unique(yo3$intid); b<-Bgroups; b<-as.character(b)
rowcount<-1
for(i in 1:length(b)){
spp<-yo3[yo3$intid==b[i],]
spp2<-spp[1:3,]
yo<-aggregate(spp2["phenodiff"], by=spp2[c("studyid","intid")], FUN=mean)
asdf<-rowcount+(nrow(spp)-1)
best[rowcount:asdf,1]<-b[i]
best[rowcount:asdf,2]<-spp[,c("year")]
best[rowcount:asdf,3]<- with(spp, phenodiff-yo[1,3])
rowcount<-rowcount+nrow(spp)
}
yo4<-merge(yo3, best, by=c("intid","year"))
yo3<-yo4

yo3$phenofreq<-as.factor(yo3$phenofreq)
yo3$taxatype<-as.factor(yo3$taxatype)
yo3$interaction[yo3$interaction=="poll"] <- "pollination"
yo3$interaction<-as.factor(yo3$interaction)


mat<-subset(yo3, length>4 & intid!="184" & intid!="186" & intid!="187" & intid!="188")