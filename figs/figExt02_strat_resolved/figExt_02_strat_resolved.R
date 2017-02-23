rm(list=ls())

########
## An example of a time-calibrated phylogeny with randomly resolved polytomies
##
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## January 2017
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
## January 2017
########

########
#bigEye#
########

#plot for SI: one example of a fully resolved tree, time-calibrated, including the stratigraphic ranges, all plotted against the geologic time scale
#note: necessary to re-timecalibrate the tree anew because write.tree() does not preserve the root age estimate which is needed to align the tree with the geologic timescale 


#preliminaries
##make sure that all required files are in your working directory and that all packages are installed

#required files
# bigEye_topology.tre
# sampled.alt.csv


#call libraries
require(ape)
require(strap)
require(paleotree)


######### STEP 1load tree and time data
phy <- read.tree("../data/paleo/bigEye_topology.tre") #make sure to always call the most recent version!
#phy <- ladderize(phy)
timeData <- read.csv("../data/paleo/sampled.alt.csv", header=T); rownames(timeData) <- timeData$taxa #these are the stratigraphic ranges
timeData <- as.matrix(timeData[,2:3])
head(timeData) #making sure everything is correct

#time-scaling the tree but leaving the polytomies
phy <- timePaleoPhy(phy, timeData, type = "equal", vartime = 1, randres = T, ntrees=1, dateTreatment= "minMax")
phy <- rotate(phy, 80)
#making sure the root time is there as it's required for alignment with geologic time scale! (lost when saved to file in tre-format...)
phy$root.time #million years

pdf("figExt_02_strat_resolved.pdf", 8, 11)
geoscalePhylo(tree=phy, ages=timeData, boxes="Age", cex.tip=0.7, edge.color="black", cex.ts=0.8, cex.age=0.8)
dev.off()

