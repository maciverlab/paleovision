rm(list=ls())

########
## Summary of our phylogenetic comparative study
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

#call libraries
require(ape)
require(geiger)
require(nlme)
require(phytools)
require(strap)
require(paleotree)

#designing figure for bayOU and relative eye size

#need:bigEye_topology.R; sampled.alt.csv; BMres.csv

#load tree and data
phy <- read.tree("../data/paleo/bigEye_topology.tre") #make sure to always call the most recent version!
timeData <- read.csv("../data/paleo/sampled.alt.csv", header=T); rownames(timeData) <- timeData$taxa #these are the stratigraphic ranges
timeData <- as.matrix(timeData[,2:3])

#eye data
eyes <- read.csv("../data/paleo/BMres.csv", header=T)
rownames(eyes) <- phy$tip.label
eyes <- eyes[phy$tip.label,]

#time-scaling the tree
phy1 <- timePaleoPhy(phy, timeData, type = "equal", vartime = 1, randres = F, ntrees=1, dateTreatment= "minMax")
phy1 <- rotate(phy1, 80)
eyes <- eyes[phy1$tip.label,]
phy1 <- reorder(phy1, order = "postorder")

#checking what branches feature shifts in fully resolved tree: from bayOU: 61, 78, 80, 92
phy2 <- timePaleoPhy(phy, timeData, type = "equal", vartime = 1, randres = T, ntrees=1, dateTreatment= "minMax")
phy2 <- reorder(phy2, order = "postorder")
#plot(phy2, cex=0.5); edgelabels(cex=0.5)

#finding edges in tree with polytomies that correspond to the bayou findings
#phy1 <- reorder(phy1, order = "postorder")
#plot(phy1, cex=0.5); edgelabels(cex=0.5)

#58, 75, 77, 85


colr <- rep("black", dim(phy$edge)[1]) 
colr[58] <- "blue"
colr[75] <- "red"
colr[77] <- "red4"
colr[85] <- "red4"

thick <- rep(1, dim(phy$edge)[1]) 
thick[c(58, 77, 85)] <- 3
thick[75] <- 5

res <- eyes[,2]*10



############################
pdf("fig02.pdf", 8, 12, useDingbats=F)
geoscalePhylo(tree=phy1, boxes="Epoch", units=c("Period", "Epoch"), label.offset=3.75, tick.scale=10,
              cex.tip=0.6, edge.color=colr, width=thick, cex.ts=1, cex.age=1)

dat <- res

circle.size <- (dat+abs(min(dat))+1)/2.5

tiplabels(pch=21, bg="dark grey", cex=circle.size, adj=2)

legend(1,60, # places a legend at the appropriate place 
       c("1.53x", "1.42x", "1.06x", "no change supported"), # puts text in the legend 
       bty="n",
       lty=c(1,1,1,1), # gives the legend appropriate symbols (lines)
       
       lwd=c(5,3,3,1),col=c("red", "red4", "blue", "black"),
       
       title="optimal eye size compared to root" ) 

#rect(2, 38, 26, 42, border="black")
text(1, 51, label="residual eye socket size", pos=4)
points(6, 50, pch=21, col="black", bg="dark grey", cex=min(circle.size))
points(6, 49, pch=21, col="black", bg="dark grey", cex=max(circle.size))
text(9.5, 49.75, label="min", pos=4)
text(9.5, 48.75, label="max", pos=4)

dev.off()
#############################
















