rm(list=ls())

########
## Corresponding adaptive peaks in relative eye sizes
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz 
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## January 2017
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
########

#summarizing information from the magnitude of theta2 along branches (bayOU results)

#required files:
# all PP/theta output from bayOU, set up to be saved in the following format: bayOU_[<tree number>]_[chain <a> or <b>]_posteriors


##make sure that all required files are in your working directory and that the MASS package is installed


require(MASS)
#read in the PP

readPP <- function (file) {
  return(tryCatch(read.csv(file), error=function(e) NULL))
}

dataFiles <- lapply(Sys.glob("../data/paleo/bayou_output/*posteriors.csv"), readPP)

dataFiles <- dataFiles[!sapply(dataFiles,is.null)]

length(dataFiles)

#make an empty matrix filled with NA

all.mag <- matrix(NA, nrow = 116, ncol = 1668) #800(?) chains were returned in total, and we have 116 branches

#read in the magnitude of theta2

for (i in 1:1668) {
  temp <- dataFiles[[i]]$magnitude.of.theta2
  all.mag[,i] <- temp
}




#calculate mean for each branch (row)
mag <- apply(all.mag, 1, mean)

#check the interesting branches (see summary of posterior probabilities)
#61 78 80 92
mag[61]
mag[78]
mag[80]
mag[92]


#check if there are any that would predict a different direction!
which(all.mag[61,]<0) 
which(all.mag[78,]<0)
which(all.mag[80,]<0)
which(all.mag[92,]<0)


pdf("figExt_03f.pdf", 4, 4)
hcols <- rgb(0, 0, 255, 150, maxColorValue=255)
hist(all.mag[61,], col=hcols, xlab="New Optimum", ylab="Density", xlim=c(0,3), ylim=c(0,500), main="")
#abline(v=mag[61], lty=2, col="red")

hcols <- rgb(255, 0, 0, 150, maxColorValue=255)
hist(all.mag[78,], col=hcols, xlab="new theta", ylab="density", add=T)
#abline(v=mag[78], lty=2, col="red")

hcols <- rgb(139, 0, 0, 150, maxColorValue=255)
hist(all.mag[80,], col=hcols, xlab="new theta", ylab="density", add=T, breaks=20)
#abline(v=mag[80], lty=2, col="red")

hcols <- rgb(139, 0, 0, 75, maxColorValue=255)
hist(all.mag[92,], col=hcols, xlab="new theta", ylab="density", add=T)

abline(v=0.2058466, lty=2, col="red")
dev.off()

###############



