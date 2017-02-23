rm(list=ls())

########
## Gelman’s R plots for log likelihood, alpha , and variance estimates from one example tree in our sample
## Title: A massive increase in visual range preceded the origin of terrestrial vertebrates
## Author: Lars Schmitz 
## Affiliation: Claremont McKenna, Pitzer, and Scripps Colleges
## January 2017
## This work is licensed under the Creative Commons Attribution 4.0 International License.
## To view a copy of this license, visit http://creativecommons.org/licenses/by/4.0/.
########


#calling libraries
require(ape)
require(phytools)
require(geiger)
require(denstrip)
require(ouch)
require(bayou)
require(MASS)
require(nlme)

##############loading data and trees
#loading data
eyes <- as.data.frame(read.csv("../data/paleo/bigeye0527.csv", header=T))
rownames(eyes) <- eyes$labels #naming each row using taxon names (required for comparing phylogeny and data; see below)
dat <- eyes

#preparing data
dat[, 2:5] <- lapply(dat[, 2:5], as.character)
dat[, 2:4] <- lapply(dat[, 2:4], as.numeric) #if there is a more elegant way to make sure the columns in the dataframe are recognized as numbers let me know...
dat[, 2:3] <- lapply(dat[, 2:3], log10)

#loading the entire tree set
tree.set <- read.tree("../data/paleo/alt.equal.trees.tre")

#selecting a tree sample (again let's keep this managable for this test run)
tree <- tree.set[[496]]

dat <- dat[tree$tip.label,]

trait <- dat

#calculating bm residuals
bmregr <- gls(om~ppl,data=trait,correlation=corPagel(0.5,tree,fixed=FALSE),method="ML")
res <- as.vector(bmregr$residuals)
names(res) <- tree$tip.label
res <- res*10 #the residuals are very small, especially given the time span of our tree (>120ma); multiplying by 10 makes likelihood calculations more stable (not running up against infinite)
#this was suggested to me by Josef Uyeda (author of the bayOU paper)


#prior
prior8 <- make.prior(tree, dists=list(dalpha="dlnorm", 
                                      dsig2="dlnorm",
                                      dsb="dsb", 
                                      dk="cdpois", 
                                      dtheta="dnorm",
                                      dloc="dunif"),
                     param=list(dk=list(lambda=12, kmax=116),dtheta=list(mean=mean(res), sd=4),dalpha=list(meanlog=0, sdlog=1), dsig2=list(meanlog=0, sdlog=1))
                     
)

fit1 <- bayou.mcmc(tree, res, SE=0.001, model="OU", prior8, ngen=400000, new.dir=getwd(), plot.freq=NULL, ticker.freq=10000)
chain <- load.bayou(fit1, save.Rdata=T, cleanup=T)
chain <- set.burnin(chain, 0.3)
out <- summary(chain)

fit2 <- bayou.mcmc(tree, res, SE=0.001, model="OU", prior8, ngen=400000, new.dir=getwd(), plot.freq=NULL, ticker.freq=10000)
chain2 <- load.bayou(fit2, save.Rdata=T, cleanup=T)
chain2 <- set.burnin(chain2, 0.3)
out2 <- summary(chain2)

L1 <- Lposterior(chain,tree, burnin=0.3)
L2 <- Lposterior(chain2,tree, burnin=0.3)

#phenogram.density(tree, dat, chain=chain2, burnin=0.3, pp.cutoff=0.3)
pdf("figExt_03abc.pdf")
par(mfrow=c(2,2))
RlnL <- gelman.R("lnL", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=120000, lty=2, col="red")
Ralpha <- gelman.R("alpha", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=120000, lty=2, col="red")
Rsig2 <- gelman.R("sig2", chain1=chain, chain2=chain2, plot=TRUE, type="n", ylim=c(0.9, 2))
abline(v=120000, lty=2, col="red")
plot(L1$pp,L2$pp, xlim=c(0,1), ylim=c(0,1), xlab="Chain 1", ylab="Chain 2")
curve(1*x, add=TRUE, lty=2)
dev.off()





