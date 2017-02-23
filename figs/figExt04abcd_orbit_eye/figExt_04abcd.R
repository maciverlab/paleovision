rm(list=ls())

########
## Eye socket eye correlation
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz 
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## January 2017
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
########

#exploring the correlation between the bony eyesocket and eye diameter in teleost fish

#required files:
# mean_ol_ed_pd_0502.csv
# fish.tre #source: Rabosky et al 2013 Nature Communications


##make sure that all required files are in your working directory and that the MASS package is installedrm(list=ls())


###load packages for R
require(ape)
require(geiger)
require(nlme)
require(MASS)

###data
eyedata <- as.data.frame(read.csv("../data/paleo/eye_orbit.csv", header=T))
ratio <- (eyedata$ed/eyedata$ol)*100
eyedata[2:3] <- lapply(eyedata[2:3],log10)


###standard (non-phylogenetic) statistics
fit.all <- lm(ed~ol, eyedata)
fit.all.sum <- summary(fit.all)

fit.all.sum$r.squared #r squared
fit.all.sum$coefficients[[2]] #slope estimate
fit.all.sum$coefficients[[4]] #se of slope

range(ratio)
mean(ratio)
var(ratio)
sd(ratio)


###phylogenetic method

#tree
tree <- read.tree("../data/paleo/fish.tre") #source: Rabosky et al 2013 Nature Communications


#trimming the tree, and preparing data
rownames(eyedata) <- eyedata$taxon
compare <- treedata(tree, eyedata, sort=T, warnings=T) # the treedata function compares data with tree and prunes the tree automatically
tree <- compare$phy
tree <- multi2di(tree)
trait <- as.data.frame(compare$data)
trait[, 2:3] <- lapply(trait[, 2:3], as.character)
trait[, 2:3] <- lapply(trait[, 2:3], as.numeric)
variables <- trait[,2:3]


#prepare data
eyecorr<-as.data.frame(cbind(taxon=tree$tip.label,variables))


#regression with phylo
fit <- lm(ed~ol, data=eyecorr) #standard least square regression for comparison
bmfit <- gls(ed~ol,data=eyecorr,correlation=corPagel(0.5,tree,fixed=FALSE),method="ML") 
#this may give warnings of false convergence. If that's the case modify the starting parameter value of lambda
#do the same if encountering a NA.NaN/Inf error
oufit <- gls(ed~ol,data=eyecorr,correlation=corMartins(0.5,tree,fixed=FALSE),method="ML")

fit.sum <- summary(fit)
bmfit.sum <- summary(bmfit)
oufit.sum <- summary(oufit)

fit.sum$r.squared #r squared
fit.sum$coefficients[[2]] #slope estimate
fit.sum$coefficients[[4]] #se of slope

bmfit.sum$tTable[[2]] #slope estimate
bmfit.sum$tTable[[4]] #se of slope

oufit.sum$tTable[[2]] #slope estimate
oufit.sum$tTable[[4]] #se of slope



###plotting combined
#1-3
pdf("figExt_04abcd.pdf", useDingbats = F)
par(mfrow=c(2,2))
plot(ed~ol, eyedata, xlab="log10 Orbit Length, mm", ylab="log10 Eye Diameter, mm", main="Eye vs. Orbit")
abline(fit.all)

boxplot(ratio, col="steel blue", horizontal=T, xlab="Eyeball Diameter/Orbit Length, %", main="Boxplot of Eye-Orbit-Ratio")

plot(ed~ol, data=eyecorr, xlab="log10 Orbit Length, mm", ylab="log10 Eye Diameter, mm", main="Eye vs. Orbit")
abline(fit)
abline(bmfit, lty=2, col="red")
abline(oufit, lty=2, col="blue")

hist(ratio, xlab="Eyeball Diameter/Orbit Length, %", col="steel blue", main="Histogram of Eye-Orbit-Ratio")

dev.off()

######

