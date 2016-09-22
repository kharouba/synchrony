#load data
rm(list=ls())
setwd("/users/kharouba/google drive/UBC/synchrony project/analysis")

mom<-read.csv("tempsens_nov3.csv", header=TRUE, na.strings="NA", as.is=TRUE)
mom<-na.omit(mom)
mom$spp<-as.factor(mom$spp)
mom2<-subset(mom, studyid!="HMK033" & studyid!="HMK030" & studyid!="SET005" & studyid!="HMK035" & intid!="182" & studyid!="EMW001" & intid!="184" & intid!="186" & intid!="187" & intid!="188")
mom<-mom2
mom2<-unique(mom[,c("studyid","year","envvalue","z","species","phenovalue","slope","intid","int_type","spp")])
mom<-mom2


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

#load yo3 from stats_2015.R
yumm<-merge(mom, unique(yo3[,c("studyid","intid","year")]), by=c("studyid","intid","year"))
mom<-yumm

ggplot(data=sub, aes(y=phenovalue, x=z, group=1))+
geom_line()+geom_point()+theme_bw()

doPlot <- function(sel_name) {
   subby <- mom[mom$species == sel_name,]
   ggobj <- ggplot(data=subby, aes(z, phenovalue)) +
       geom_point(size=2) + 
       #geom_line()+
       geom_smooth(method="lm", se=FALSE, aes(colour="red"))+stat_smooth(method = "lm", formula = y ~ x + I(x^2), colour="blue", se=FALSE)+
       theme_bw()+theme(legend.position="none",axis.title.x =element_text(size=17), axis.text.x=element_text(size=17), axis.text.y=element_text(size=17), axis.title.y=element_text(size=17, angle=90))+ylab("Day of year")+xlab("z")
   print(ggobj)
   ggsave(sprintf("/users/kharouba/google drive/UBC/synchrony project/graphs/doybyz/species/spp%spp.pdf", sel_name))
}
   
lapply(unique(mom$species), doPlot)


plot(phenovalue~z, data=sub, type="l")

        points(phenovalue~year, data=subby.nohinge, cex=0.6)

        abline((lm(phenovalue~year, data=subby.nohinge)))

        lines(phenovalue~newyear, data=subby.hinge, col="red")

        abline((lm(phenovalue~newyear, data=subby.hinge)), col="red")
        
        
        
  
