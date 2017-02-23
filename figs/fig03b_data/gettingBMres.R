rm(list=ls())

########
## Bayou residulas for eye size
##
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
## January 2017
########

########
#bigEye#
########
########
#bigEye#
########


#call libraries
require(ape)
require(geiger)
require(nlme)
require(phytools)


#load tree and data
phy.all <- read.tree("../data/paleo/alt.equal.trees.tre") #reading in tree samples
eyes <- read.csv("../data/paleo/bigeye0527.csv", header=T)


#set up empty vectors for storing results

all.res <- matrix(NA, nrow = 59, ncol = 1000)


##############################start the loop
for (i in 1:1000) try({


###########preliminaries
phy <- phy.all[[i]] #just go with one single tree for now
#phy2 <- phy
#phy2$edge.length <- phy2$edge.length*10
#phy$edge.length <- phy$edge.length*10
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

BM.scaling <- gls(om~ppl,data=bigeye,correlation=corPagel(0.5,phy,fixed=FALSE),method="ML")

all.res[,i] <- BM.scaling$residuals
#
})



##############end of loop


res <- apply(all.res, 1, mean)
names(res) <- phy.all[[1]]$tip.label
write.csv(res, "BMres.csv")
