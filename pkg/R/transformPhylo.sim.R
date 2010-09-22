transformPhylo.sim <- function(phy, n=1, x=NULL, model=NULL, kappa=NULL, lambda=NULL, delta=NULL, alpha=NULL, psi=NULL, nodeIDs=NULL, rateType=NULL, cladeRates=NULL, branchRates=NULL, rate=NULL, group.means=NULL) {
	
	switch(model,		  
		   
		   "bm" = {
					transformPhy <- phy
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
		   
		   "kappa" = {
					transformPhy <- transformPhylo(phy=phy, model="kappa", kappa=kappa)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
		   
		   "lambda" = {
					transformPhy <- transformPhylo(phy=phy, model="lambda", lambda=lambda)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
		   
		   "delta" = {
					transformPhy <- transformPhylo(phy=phy, model="delta", delta=delta)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
					
		   "free" = {
					transformPhy <- transformPhylo(phy=phy, model="free", branchRates=branchRates)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)		
					},
		   
		   "clade" = {
					transformPhy <- transformPhylo(phy=phy, model="clade", nodeIDs=nodeIDs, cladeRates=cladeRates, rateType=rateType)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
		   
		   "OU" = {
					transformPhy <- transformPhylo(phy=phy, model="OU", alpha=alpha)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
		   "psi" = {
					transformPhy <- transformPhylo(phy=phy, model="psi", psi=psi)
					phyMat <- vcv.array(transformPhy)
					attr(phyMat, "class") <- "matrix"
					ydum <- as.matrix(t(rmvnorm(n, sigma = phyMat)))
					rownames(ydum) <- rownames(phyMat)
					},
		   
		   "mixedRate" = {
					
					x <- as.matrix(x)
					dat <- data.frame(x=x, y=rep(0, length(x[,1])))
					
					rateData <- as.rateData(y=dat$y, x=dat$x, rateMatrix = NULL, phy=phy, data=dat)
		   
					V <- transformRateMatrix(rateData, rate=rate)
					expect.sd <- sqrt(mean(V[upper.tri(V)]))
		   
					if (is.null(group.means)) {ydum <- as.matrix(t(rmvnorm(n, sigma =  (V) )))
						rownames(ydum) <- rownames(V)} 
		   
					else {
						x.means <- unique(rateData$x)
						n.means <- length(x.means)
						samp.means <- rep(NA, length(rateData$x))
						ydum <- vector(mode="list", length=length(group.means))
						for (i in 1:n.means) {
							samp.means[which(rateData$x == (i-1))] <- rep(0+(expect.sd*group.means[i]), length(which(rateData$x == (i-1))))
							}	
		   
						ydum <- as.matrix(t(rmvnorm(n, mean=samp.means, sigma =  (V) )))
						rownames(ydum) <- rownames(V)
						}
					}
		   )
	
	
	return(ydum)
}
			