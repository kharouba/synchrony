# set up as much as you can in advance that will be in the figure, that means, data
to use, models to use etc., definitely set up vectors to hold colours and pch to
step through

man<-read.csv("rough_man.csv", header=TRUE, na.strings="NA", as.is=TRUE)
daterr <- man
siteIDname <- unique(daterr$studyid)
colourz <- with(daterr, rainbow(length(unique(studyid))))
pchobs <- c(rep(1:11,4))


quartz("Quartz", width=8, height=5, pointsize=12) # open a figure

# make a plot without any points (but getting the x and y axes large enough for all
the data)
with(daterr, plot(fitnessvalue~treatchangeday, type="n",
    ylab="fitness", xlab="change in days"))
# loop through each study to add in points
for(i in 1:length(siteIDname)){
   dater <- subset(daterr, studyid==siteIDname[i])
   with(dater, points(fitnessvalue~treatchangeday, col=colourz[i], pch=pchobs[i]))
# can add regression and abline here, be sure to set the colour on the abline
}


## this is the code (below) behind figure 4 in my meta-analysis, we had 3 different
models to look at so the code builds on those

######
## meanevent plot, pretty
######

mod.flo.event <- lme(sens~meaneventdate*datype, random=~1|site/latbi,
data=subset(sens.m.nl,
     phentype=="flo"), na.action=na.omit)
mod.flo.eventout <- lme(sens~meaneventdate*datype, random=~1|site/latbi,
     data=subset(sens.m.nlout, phentype=="flo"), na.action=na.omit)
mod.flo.eventout.wg <- lme(sens~meaneventdate*datype, random=~1|site/latbi,
   data=subset(sens.m.nlout, phentype=="flo" & site !="gothic"
   & site!="washdc"), na.action=na.omit)

modeluse <- mod.flo.eventout
daterr <-subset(sens.m.nlout, phentype=="flo")
ablinex
<-fixed.effects(modeluse)[["(Intercept)"]][1]+fixed.effects(modeluse)[["datypeobs"]][1]
abliney <-fixed.effects(modeluse)[["meaneventdate"]][1]+
    fixed.effects(modeluse)[["meaneventdate:datypeobs"]][1]
ablinex.exp <-fixed.effects(modeluse)[["(Intercept)"]][1]
abliney.exp <-fixed.effects(modeluse)[["meaneventdate"]][1]
leg.cex <- 0.5
legend.x1 <- 20
legend.y1 <- 33
legend.x <- 260
legend.y <- -6 # all data: 63
exp <- subset(daterr, datype=="exp")
obs <- subset(daterr, datype=="obs")
expsitenames <- unique(exp[["site"]])
obssitenames <- unique(obs[["site"]])
colobs <- rep("mediumturquoise", 50)
colexp <- c(rep("firebrick4", 15), rep("violetred4", 10))
pchobs <- c(rep(1:11,2))
pchexp <- c(rep(1:24))
par(mfrow=c(1,1))
quartz("Quartz", width=8, height=5, pointsize=12)
plot(daterr$sens~daterr$meaneventdate, type="n",
    ylab="sensitivity", xlab="mean event date")
for(i in c(1:length(obssitenames))){
   dater <- subset(obs, site==obssitenames[i])
   points(dater$sens~dater$meaneventdate, col=colobs[i], pch=pchobs[i])
 }
for(i in c(1:length(expsitenames))){
   dater <- subset(exp, site==expsitenames[i])
   points(jitter(dater$sens)~dater$meaneventdate, col=colexp[i], pch=pchexp[i])
 }
abline(ablinex, abliney, col=colobs[1])
abline(ablinex.exp, abliney.exp, col=colexp[1])
legendsitesobs <- merge(as.data.frame(obssitenames), sitenums, by.x="obssitenames",
by.y="site")

legend(legend.x, legend.y, legendsitesobs$numbers, pch=pchobs,
    col=colobs, cex=leg.cex, bty="n")
sitelist <- subset(expsites, select=c("site", "location"))
legendsites <- merge(as.data.frame(expsitenames), sitenums, by.x="expsitenames",
by.y="site")
legend(legend.x1, legend.y1, legendsites$numbers, pch=pchexp,
     col=colexp, cex=leg.cex, bty="n")
