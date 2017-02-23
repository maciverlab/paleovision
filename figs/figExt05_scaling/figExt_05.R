rm(list=ls())

########
## slope estimates and the associated standard error over all 1000 trees; OU, BM
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz 
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## January 2017
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
########

#investigating whether size correction is necessary (see more detailed explanation lines 29ff.)

#required files:
# alt.equal.trees.tre
# bigeye0527.csv

#call libraries
require(ape)
require(geiger)
require(nlme)
require(phytools)

#load tree and data
phy.all <- read.tree("../data/paleo/alt.equal.trees.tre") #reading in tree samples
eyes <- read.csv("../data/paleo/bigeye0527.csv", header=T)


#there will be three main parts with IDENTICAL approaches (PGLS, BM and OU)

#1 is the ratio size dependent?
  #a) OU slope and standard error; PLOT needed at the end
  #b) p-value
  #c) AIC score

  #d) BM slope and standard error; PLOT needed at the end
  #e) p-value
  #f) AIC score

  #####g) dAIC score to see whether one or the other fits better (not crucial though especially if both agree on slope estimates)

#2 what's the slope of om~ppl?
  #a) OU slope and standard error; PLOT needed at the end
  #b) p-value
  #c) AIC score
  
  #d) BM slope and standard error; PLOT needed at the end
  #e) p-value
  #f) AIC score
  
  ######g) dAICc score to see whether one or the other fits better (not crucial though especially if both agree on slope estimates)

#3 Was size correction successful?
  #a) OU slope and standard error; PLOT needed at the end
  #b) p-value
  #c) AICc score
  
  #d) BM slope and standard error; PLOT needed at the end
  #e) p-value
  #f) AICc score
  
  ######g) dAICc score to see whether one or the other fits better (not crucial though especially if both agree on slope estimates)

#set up empty vectors for storing results

result_1a_slope <- as.vector(rep(NA, 1000))
result_1a_se <- as.vector(rep(NA, 1000))
result_1b_p <- as.vector(rep(NA, 1000))
result_1c_AIC <- as.vector(rep(NA, 1000))
result_1d_slope <- as.vector(rep(NA, 1000))
result_1d_se <- as.vector(rep(NA, 1000))
result_1e_p <- as.vector(rep(NA, 1000))
result_1f_AIC <- as.vector(rep(NA, 1000))

result_2a_slope <- as.vector(rep(NA, 1000))
result_2a_se <- as.vector(rep(NA, 1000))
result_2b_p <- as.vector(rep(NA, 1000))
result_2c_AIC <- as.vector(rep(NA, 1000))
result_2d_slope <- as.vector(rep(NA, 1000))
result_2d_se <- as.vector(rep(NA, 1000))
result_2e_p <- as.vector(rep(NA, 1000))
result_2f_AIC <- as.vector(rep(NA, 1000))

result_3a_slope <- as.vector(rep(NA, 1000))
result_3a_se <- as.vector(rep(NA, 1000))
result_3b_p <- as.vector(rep(NA, 1000))
result_3c_AIC <- as.vector(rep(NA, 1000))
result_3d_slope <- as.vector(rep(NA, 1000))
result_3d_se <- as.vector(rep(NA, 1000))
result_3e_p <- as.vector(rep(NA, 1000))
result_3f_AIC <- as.vector(rep(NA, 1000))

##############################start the loop
for (i in 1:1000) try({


###########preliminaries
phy <- phy.all[[i]] #just go with one single tree for now
phy2 <- phy
phy2$edge.length <- phy2$edge.length*10
rownames(eyes) <- eyes$labels #assigning rownames (used for matching purposes)
ordered.eyes <- eyes[phy$tip.label,]

variables <- ordered.eyes[,2:4]
ppl0 <- variables[,1]
names(ppl0) <- phy$tip.label
om0 <- variables[,2]
names(om0) <- phy$tip.label
groups <- ordered.eyes$group
names(groups) <- phy$tip.label
ppl <- signif(log10(ppl0), 3)
names(ppl) <- phy$tip.label
om <- signif(log10(om0),3)
names(om) <- phy$tip.label
ratio <- variables[,3]
names(ratio) <- phy$tip.label

#putting it all in a dataframe
bigeye<-as.data.frame(cbind(ppl,om, ppl0, om0, ratio, groups))
bigeye<-cbind(phy$tip.label,bigeye)
colnames(bigeye)<-c("taxon","ppl","om", "ppl0", "om0", "ratio", "group")






##############1 is the ratio size dependent?
#a) OU slope and standard error; PLOT needed at the end
OU.sizeeffect.ratio <- try(gls(ratio~ppl,data=bigeye,correlation=corMartins(0.03,phy,fixed=FALSE),method="ML"))
if(class(OU.sizeeffect.ratio)=="try-error"){OU.sizeeffect.ratio <- try(gls(ratio~ppl,data=bigeye,correlation=corMartins(0.01,phy,fixed=FALSE),method="ML"))}
if(class(OU.sizeeffect.ratio)=="try-error"){remove(OU.sizeeffect.ratio); OU.sizeeffect.ratio<- try(gls(ratio~ppl,data=bigeye,correlation=corMartins(0.003,phy,fixed=FALSE),method="ML"))}

OU1a <- summary(OU.sizeeffect.ratio)
result_1a_slope[[i]] <- OU1a$tTable[[2]]
result_1a_se[[i]] <- OU1a$tTable[[4]]

#b) p-value
result_1b_p[[i]] <- OU1a$tTable[[8]]

#c) AIC score
result_1c_AIC[[i]] <- OU1a$AIC

#d) BM slope and standard error; PLOT needed at the end
BM.sizeeffect.ratio <- gls(ratio~ppl,data=bigeye,correlation=corPagel(0.5,phy,fixed=FALSE),method="ML")
BM1d <- summary(BM.sizeeffect.ratio)
result_1d_slope[[i]] <- BM1d$tTable[[2]]
result_1d_se[[i]] <- BM1d$tTable[[4]]

#e) p-value
result_1e_p[[i]] <- BM1d$tTable[[8]]

#f) AIC score
result_1f_AIC[[i]] <- BM1d$AIC




##############2 what's the slope of om~ppl?
#a) OU slope and standard error; PLOT needed at the end
OU.scaling <- try(gls(om~ppl,data=bigeye,correlation=corMartins(0.03,phy,fixed=FALSE),method="ML"))
if(class(OU.scaling)=="try-error"){OU.scaling <- try(gls(om~ppl,data=bigeye,correlation=corMartins(0.01,phy,fixed=FALSE),method="ML"))}
if(class(OU.scaling)=="try-error"){remove(OU.scaling);OU.scaling <- try(gls(om~ppl,data=bigeye,correlation=corMartins(0.001,phy2,fixed=FALSE),method="ML"))}

OU2a <- summary(OU.scaling)
result_2a_slope[[i]] <- OU2a$tTable[[2]]
result_2a_se[[i]] <- OU2a$tTable[[4]]

#b) p-value
result_2b_p[[i]] <- OU2a$tTable[[8]]

#c) AIC score
result_2c_AIC[[i]] <- OU2a$AIC

#d) BM slope and standard error; PLOT needed at the end
BM.scaling <- gls(om~ppl,data=bigeye,correlation=corPagel(0.5,phy,fixed=FALSE),method="ML")
BM2d <- summary(BM.scaling)
result_2d_slope[[i]] <- BM2d$tTable[[2]]
result_2d_se[[i]] <- BM2d$tTable[[4]]

#e) p-value
result_2e_p[[i]] <- BM2d$tTable[[8]]

#f) AIC score
result_2f_AIC[[i]] <- BM2d$AIC




######################3 Was size correction successful?

#let's get the residuals from the fitted lines above
OU_res <- as.vector(OU.scaling$residuals)
names(OU_res) <- phy$tip.label

BM_res <- as.vector(BM.scaling$residuals)
names(BM_res) <- phy$tip.label

#and put them into a  dataframe with PPL
oures.dat <- as.data.frame(cbind(OU_res=OU_res, BM_res=BM_res, ppl=ppl))

#######

#a) OU slope and standard error; PLOT needed at the end
OU.test <- try(gls(OU_res~ppl,data=oures.dat,correlation=corMartins(0.03,phy,fixed=FALSE),method="ML"))
if(class(OU.test)=="try-error"){OU.test <- try(gls(OU_res~ppl,data=oures.dat,correlation=corMartins(0.01,phy,fixed=FALSE),method="ML"))}
if(class(OU.test)=="try-error"){remove(OU.test);OU.test <- try(gls(OU_res~ppl,data=bigeye,correlation=corMartins(0.001,phy,fixed=FALSE),method="ML"))}


OU3a <- summary(OU.test)
result_3a_slope[[i]] <- OU3a$tTable[[2]]
result_3a_se[[i]] <- OU3a$tTable[[4]]

#b) p-value
result_3b_p[[i]] <- OU3a$tTable[[8]]

#c) AIC score
result_3c_AIC[[i]] <- OU3a$AIC

#d) BM slope and standard error; PLOT needed at the end
BM.test <- gls(BM_res~ppl,data=oures.dat,correlation=corPagel(0.5,phy,fixed=FALSE),method="ML")
BM3d <- summary(BM.test)
result_3d_slope[[i]] <- BM3d$tTable[[2]]
result_3d_se[[i]] <- BM3d$tTable[[4]]

#e) p-value
result_3e_p[[i]] <- BM3d$tTable[[8]]

#f) AICc score
result_3f_AIC[[i]] <- BM3d$AIC


})



##############end of loop

#save everything, inclduing junk, for complete record
results <- as.data.frame(cbind(result_1a_slope, 
                              result_1a_se ,
                              result_1b_p ,
                              result_1c_AIC, 
                              result_1d_slope, 
                              result_1d_se ,
                              result_1e_p ,
                              result_1f_AIC, 
                              
                              result_2a_slope, 
                              result_2a_se ,
                              result_2b_p ,
                              result_2c_AIC, 
                              result_2d_slope, 
                              result_2d_se ,
                              result_2e_p ,
                              result_2f_AIC, 
                              
                              result_3a_slope, 
                              result_3a_se ,
                              result_3b_p ,
                              result_3c_AIC, 
                              result_3d_slope, 
                              result_3d_se ,
                              result_3e_p ,
                              result_3f_AIC))

write.csv(results, "scaling.results.csv")


######### and now let's plot stuff


#clean up results!
junk1 <- which(result_1c_AIC>-160) #all these seem to have flat likelihood curves!!! must exclude
result_1a_slope[junk1] <- NA

junk2 <- which(result_2c_AIC>-60) #somewhere between 45 and 50 neg
result_2a_slope[junk2] <- NA

junk3 <- which(result_3c_AIC>-52) #somewhere between 45 and 50 neg
result_3a_slope[junk3] <- NA


mean(result_1a_slope, na.rm=T)
mean(result_1d_slope, na.rm=T)
mean(result_2a_slope, na.rm=T)
mean(result_2d_slope, na.rm=T)
mean(result_3a_slope, na.rm=T)
mean(result_3d_slope, na.rm=T)

#function to visualize the slopes and their standard errors over the entire tree sample
slopesum <- function (slopes, se, minslope, maxslope, title, ...){
  llim <- slopes-se
  ulim <- slopes+se
  
  plot(slopes, as.vector(1:length(slopes)), type="n", main=title, xlab="slope estimates", ylab="tree sample #", xlim=(c(minslope, maxslope)))
  
  for (i in 1:length(slopes)){
    lines(x=c(llim[[i]], ulim[[i]]), y=c(i,i), col="light grey")
  }
  
  points(slopes, as.vector(1:length(slopes)), pch=4, cex=0.1)
  
}

pdf("figExt_05_scaling.pdf", useDingbats = FALSE)

par(mfrow=c(3,2))

#result_2a_slope, result_2a_se 
slopesum(result_2a_slope, result_2a_se, 0.5, 1.2, title="OM versus PPL (OU)")
abline(v=1, lty=2, col="red")

#result_2d_slope, result_2d_se
slopesum(result_2d_slope, result_2d_se, 0.5, 1.2, title="OM versus PPL (BM)")
abline(v=1, lty=2, col="red")

#result_1a_slope, ... _se
slopesum(result_1a_slope, result_1a_se, -0.3, 0.1, title="OM/PPL versus PPL (OU)")
abline(v=0, lty=2, col="red")

#result_1d_slope, ... _se
slopesum(result_1d_slope, result_1d_se, -0.3, 0.1, title="OM/PPL versus PPL (BM)")
abline(v=0, lty=2, col="red")

#result_3a_slope, result_3a_se 
slopesum(result_3a_slope, result_3a_se, -0.2, 0.2, title="Residuals versus PPL (OU)")
abline(v=0, lty=2, col="red")

#result_3d_slope, result_3d_se 
slopesum(result_3d_slope, result_3d_se, -0.2, 0.2, title="Residuals versus PPL (BM)")
abline(v=0, lty=2, col="red")

dev.off()
