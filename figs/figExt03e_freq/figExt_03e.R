rm(list=ls())

########
## Posterior probabilities of shifts
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz 
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## January 2017
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
########
#summarizing information from the posterior probabilities for shifts along branches (bayOU results)

#required files:
# all PP output from bayOU, set up to be saved in the following format: bayOU_[<tree number>]_[chain <a> or <b>]_posteriors


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

all.pp <- matrix(NA, nrow = 116, ncol = 1668) #800(?) chains were returned in total, and we have 116 branches

for (i in 1:1668) {
  temp <- dataFiles[[i]]$pp
  all.pp[,i] <- temp
}

#calculate mean for each branch (row)

pp <- apply(all.pp, 1, mean)
which(pp>0.12)#these seem to be the branches that stand out judged from the PP distribution (see below)

pp[which(pp>0.12)]

#61 78 80 92

pdf("figExt_03e.pdf", 4, 4)
den <- density(pp, adjust=5); plot(den, xlab="Posterior Probability", main="")
cols <- rgb(200, 0,0, 150, maxColorValue=255)
polygon(den, col=cols)
abline(v=1/116, col="black", lty=2)
arrows(pp[61], 15, pp[61], 3, col="steel blue", lwd=2, length=.06)
arrows(pp[78], 15, pp[78], 3, col="steel blue", lwd=2, length=.06)
arrows(pp[80], 15, pp[80], 3, col="steel blue", lwd=2, length=.06)
arrows(pp[92], 15, pp[92], 3, col="steel blue", lwd=2, length=.06)
#abline(v=pp[61], col="blue")
#abline(v=pp[78], col="blue")
#abline(v=pp[80], col="blue")
#abline(v=pp[92], col="blue")
dev.off()
