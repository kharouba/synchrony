### Figures ###



Appendix
# Temp change figure
asdf<-summary(temp.model, pars=c("a","b"))
lab<-as.data.frame(asdf[[1]][1:18]); names(lab)[1]<-"int"
lab$slope<-asdf[[1]][19:36]
lab$id<-unique(species)

me<-summary(temp.model, pars=c("a_study"))
mean(me[[1]][1:13]) #rough guess on mu_a

temp.id3$year2<-with(temp.id3, year-1981)
temp.id3$lowerCI<-0.08*temp.id3$year2+(1.26-0.06)
temp.id3$upperCI<-0.08*temp.id3$year2+(1.26+0.1)

b<-ggplot(temp.id3, aes(y=envvalue,x=year2, colour=factor(datasetid)))+geom_line()+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab(expression(paste("Temperature ",degree,"C")))+theme(legend.position="none")+xlab("Year")+geom_vline(xintercept=0, linetype="dashed", size=0.3)+geom_abline(slope=0.08, intercept=1.26, size=0.75)+geom_ribbon(aes(x=year2, ymin=lowerCI, ymax=upperCI), colour="grey", fill="grey", alpha=0.6)+annotation_custom(grob = textGrob(label = "b)", hjust = 0, gp = gpar(cex = 1.5)), ymin = 19, ymax = 19, xmin = -13, xmax = -13)

+geom_abline(aes(intercept=int, slope=slope, colour=factor(id)), data=lab)



#temp sens figure
asdf<-summary(tempsens.model, pars=c("a","b"))
lab<-as.data.frame(asdf[[1]][1:37]); names(lab)[1]<-"int"
lab$slope<-asdf[[1]][38:74]
lab$id<-unique(species)
clim3$ymin<--6.65
clim3$ymax<--3.02
clim3$lowerCI<--4.79*clim3$envvalue+(147-6.65)
clim3$upperCI<--4.79*clim3$envvalue+(147+3.00)

me<-summary(tempsens.model, pars=c("a"))
mean(me[[1]][1:37]) #rough guess on mu_a

suck<-summary(tempsens.model, pars="ypred")
nope<-as.data.frame(suck[[1]][1:764]); names(nope)[1]<-"ypred"
#nope$grp<-rep(1:37, each=1000); names(nope)[1]<-"ypred"
#yes <- aggregate(nope["ypred"], nope["grp"], FUN=mean); 
#data2<-merge(yes, data, by="grp")


a<-ggplot(clim3, aes(y=phenovalue,x=envvalue, colour=factor(species)))+geom_point()+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("Phenology (doy)")+theme(legend.position="none")+xlab(expression(paste("Temperature ",degree,"C")))+geom_abline(aes(intercept=int, slope=slope, colour=factor(id)), data=lab)+geom_abline(slope=-4.78, intercept=147, size=1)+geom_ribbon(aes(x=envvalue, ymin=lowerCI, ymax=upperCI), colour="grey", fill="grey", alpha=0.6)+annotation_custom(grob = textGrob(label = "a)", hjust = 0, gp = gpar(cex = 1.5)), ymin = 335, ymax = 335, xmin = -8.5, xmax = -8.5)

multiplot(a,b, cols=2)


# Hinge vs. non-hinge

sub<-subset(rawlong.tot, intid==170)
sub$int<-c(rep(138.07, 37), rep(98.19, 37))
sub$slope<-c(rep(-0.80, 37), rep(-0.62, 37))
#1981 hinge:
Daphnia3 spp. (spp 29); a=138.07; b=-0.80
Phytoplankton1 spp. (66); a=98.19 ; b=-0.62
# no- hinge
Daphnia3 spp. (spp 29); a=138.05 (95%: 132.85,143.24); b= -0.80 (-1.21,-0.39)
Phytoplankton1 spp. (66); a=98.18 (92.91,103.39); b-0.52 (-0.83, -0.20)

#2nd example
sub<-subset(rawlong.tot, intid==171)
sub$year2<-with(sub, year-1981)
#1981 hinge:
Daphnia3 spp. (spp 29); a=138.07 (132.90, 143.24); b=-0.80 (-1.20,-0.39)
Perca fluviatillis (64); a=139.71 (134.49, 144.86) b=-0.55 (-0.96,-0.15)
#nohinge
Daphnia3 spp. (spp 29); a=135.43 (); b= -0.66 ()
Perca fluviatillis (64); a=137.90 () b=-0.64 ()

TRY GEOM_LINE
WITH GEOM_RIBBON(aes(ymin=, ymax=), alpha=0.2)


a<-ggplot(sub, aes(y=phenovalue,x=year2, colour=factor(species)))+geom_point()+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("Phenology (doy)")+theme(legend.position="none")+xlab("")+scale_size_continuous(guide=FALSE)+geom_line(aes(colour=factor(species)))+geom_segment(aes(x=-12, xend=0, y=138.07+0, yend=138.07+0*0, linetype='hinge'),colour='#F8766D')+geom_segment(aes(x=0, xend=27, y=139.07+-(0.80), yend=138.07+(-0.80*27), linetype='hinge'),colour='#F8766D')+geom_abline(aes(slope=-0.66, intercept=135.43,linetype='no hinge'),colour='#F8766D')+geom_segment(aes(x=-12, xend=0, y=139.6+0, yend=139.6+0*0, linetype='hinge'),colour='#00BFC4')+geom_segment(aes(x=0, xend=27, y=139.71+(-0.55), yend=139.71+(-0.55*27), linetype='hinge'),colour='#00BFC4')+geom_abline(aes(slope=-0.64, intercept=137.90,linetype='no hinge'), colour='#00BFC4') +scale_linetype_manual(values=c('solid','dashed'))+xlim(-20,27)+geom_vline(xintercept=0, linetype="dotted")
# to check colours: ggplot_build(a)$data
#geom_ribbon(aes(ymin=, ymax=), alpha=0.2)

#3rd example
sub2<-subset(rawlong.tot, intid==180)
sub2$year2<-with(sub2, year-1981)
#1981 hinge:
Operophtera brumata (52) a=160.21, b=-0.61
Parus2 major (58) a=122.93 b=-0.43
#no hinge:
Operophtera brumata (52): a=156.04; b=-0.35
Parus2 major (58): a=119.92; b=-0.24

#hinge
b<-ggplot(sub2, aes(y=phenovalue,x=year2))+geom_point(aes(colour=factor(species)))+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+theme_bw()+ylab("Phenology (doy)")+xlab("Year")+geom_line(aes(colour=factor(species)))+scale_size_continuous(guide=FALSE)+geom_segment(aes(x=-20, xend=0, y=160.21+0, yend=160.21+0*0, linetype='hinge'),colour='#F8766D')+geom_segment(aes(x=0, xend=25, y=160.21+(-0.61), yend=160.21+(-0.61*25), linetype='hinge'),colour='#F8766D')+geom_abline(aes(slope=-0.35, intercept=156.04,linetype='no hinge'),colour='#F8766D')+geom_segment(aes(x=-20, xend=0, y=122.93+0, yend=122.93+0*0, linetype='hinge'),colour='#00BFC4')+geom_segment(aes(x=0, xend=25, y=122.93+(-0.43), yend=122.93+(-0.43*25), linetype='hinge'),colour='#00BFC4')+geom_abline(aes(slope=-0.24, intercept=119.92,linetype='no hinge'),colour='#00BFC4')+scale_linetype_manual(values=c('solid','dashed'))+xlim(-20,27)+geom_vline(xintercept=0, linetype="dotted")+theme(legend.position="none")

multiplot(a,b,cols=1)

# MAP  #
library(rworldmap)
newmap<-getMap(resolution="coarse")
dad<-read.csv("/users/kharouba/google drive/UBC/synchrony project/analysis/stan_2016/input/studies_geog.csv", header=TRUE, na.strings="NA", as.is=TRUE)
coordinates(dad)<-c("x","y")
plot(newmap)
plot(dad, add=T, pch=20, col="red")


#Repeating species
uni<-unique(rawlong.tot2[,c("studyid","species")])
uni$count<-1
uni_study <- aggregate(uni["count"],uni[c("studyid")], FUN=length)

uni<-unique(rawlong.tot2[,c("species","studyid")])
uni$count<-1
uni_spp<- aggregate(uni["count"],uni[c("species")], FUN=length)



within<-ggplot(uni_study, aes(x=count))+geom_histogram(binwidth=.5, alpha=.5, position="identity", colour="black")+theme(axis.title.x = element_text(size=15), axis.text.x=element_text(size=15), axis.text.y=element_text(size=15), axis.title.y=element_text(size=15, angle=90))+ylab("Number of species per study")+xlab("Number")+theme_bw()



so<-subset(rawlong.tot, spp=="spp1"); spp1<-unique(so[,c("studyid","intid","figcode")]); names(spp1)[3]<-"figcode_1"
sp<-subset(rawlong.tot, spp=="spp2"); spp2<-unique(sp[,c("studyid","intid","figcode")]); names(spp2)[3]<-"figcode_2"
tots<-merge(spp1,spp2, by=c("studyid","intid"))
comp<-merge(andtheanswer, tots, by=c("studyid","intid"))

re_uni<-ggplot(comp, aes(x=model.x, fill=figcode_1)) + 
     geom_histogram(data =subset(comp, figcode_1=="unique") , fill = "grey", col="black", alpha=0.2)+theme_bw()+ggtitle("Resource_unique")+xlim(-1.5,1.5)+ylim(0,5)
re<-ggplot(comp, aes(x=model.x, fill=figcode_1))+
     geom_histogram(data = subset(comp, figcode_1!="unique"))+theme_bw()+ggtitle("Resource_repeat")+theme(legend.position="none")+xlim(-1.5,1.5)+ylim(0,5)
co_uni<-ggplot(comp, aes(x=model.x, fill=figcode_2)) + 
     geom_histogram(data =subset(comp, figcode_2=="unique") , fill = "grey", col="black", alpha=0.2)+theme_bw()+ggtitle("Consumer_unique")+theme(legend.position="none")+xlim(-1.5,1.5)
co<-ggplot(comp, aes(x=model.x, fill=figcode_2)) + 
     geom_histogram(data = subset(comp, figcode_2!="unique"))+theme_bw()+ggtitle("Consumer_repeat")+xlim(-1.5,1.5)+theme(legend.position="none")+ylim(0,6)
multiplot(re_uni, co_uni, re, co, cols=2)
