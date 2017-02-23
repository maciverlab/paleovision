rm(list=ls())

########
## MCMC chains
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


# read in the PP

read <- function (file) {
  return(tryCatch(read.csv(file), error=function(e) NULL))
}

dataFiles1 <- lapply(Sys.glob("../data/paleo/bayou_output/*chain1_posteriors.csv"), read)

dataFiles1 <- dataFiles1[!sapply(dataFiles1,is.null)]

length(dataFiles1)

results1 <- list()

chains1.list <- for (i in 1:length(dataFiles1)) {
  temp <- dataFiles1[[i]]$pp
  results1[[i]] <- temp
  }

chains1 <- unlist(results1)

############################################
dataFiles2 <- lapply(Sys.glob("../data/paleo/bayou_output/*chain2_posteriors.csv"), read)

dataFiles2 <- dataFiles2[!sapply(dataFiles2,is.null)]

results2 <- list()

chains2.list <- for (i in 1:length(dataFiles2)) {
  temp <- dataFiles2[[i]]$pp
  results2[[i]] <- temp
}

chains2 <- unlist(results2)


###########################################

#plot(chains1, chains2)
all.colors <- colors()
selection <- rev(all.colors[173:253])

pdf("figExt_03d.pdf", useDingbats=F, 4, 4)
smoothScatter(chains1, chains2, nbin=100, nrpoints=50, xlab="Posterior Probability from Chain 1", ylab="Posterior Probability from Chain 2", bandwidth=0.005, colramp = colorRampPalette(c("white", selection)))
#plot(chains1, chains2, pch=20, cex=0.4, xlab="Posterior Probability from Chain 1", ylab="Posterior Probability from Chain 2")
curve(1*x, add=TRUE, col="red", lty=2)
#fit1 <- lm(chains2~chains1)
#abline(fit1)
#fit1
dev.off()
