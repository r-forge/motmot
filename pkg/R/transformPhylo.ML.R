transformPhylo.ML <- function(y, phy, model=NULL, nodeIDs=NULL, minCladeSize=10, nSplits=10, rateMin=1e-5, rateMax=1e5, deltaMax=10){
	
	switch(model,
		   
		   "bm" = {
					vo <- likTraitPhylo(y, phy)
					},
		   
		   "kappa" = {
					kappa <- 1
					var.funkappa <- function(kappa) { return(transformPhylo.ll(y, phy, kappa, model="kappa")[[2]])}
					vo <- optimise(var.funkappa, interval=c(0,1), maximum=TRUE)
					#vo <- optim(kappa, var.fun,  method="L-BFGS-B", lower=1e-5, upper=1, control=c(fnscale=-1))
					},

		   "lambda" = {
					lambda <- 1
					var.funlambda <- function(lambda) { return(transformPhylo.ll(y=y, phy=phy, lambda=lambda, model="lambda")[[2]])}
					vo <- optimise(var.funlambda, interval=c(0,1), maximum=TRUE)
					},
		   
		   "delta" = {
					delta <- 1
					var.fundelta <- function(delta) { return(transformPhylo.ll(y=y, phy=phy, delta=delta, model="delta")[[2]])}
					vo <- optimise(var.fundelta, interval=c(0,deltaMax), maximum=TRUE)
					},	
					
		   "free" = {
					branchRates <- rep(1,length(phy$edge.length))
					var.funfree <- function(branchRates) { return(transformPhylo.ll(y, phy, branchRates=branchRates, model="free")[[2]])}
					vo <- optim(branchRates, var.funfree, method="L-BFGS-B", lower=rateMin, upper=rateMax, control=c(fnscale=-1))
					},
		   
		   "clade" = {
					cladeRates <- rep(1,length(nodeIDs))
					var.funclade <- function(cladeRates) { return(transformPhylo.ll(y, phy, nodeIDs=nodeIDs, cladeRates=cladeRates, model="clade")[[2]])}
					vo <- optim(cladeRates, var.funclade, method="L-BFGS-B", lower=rateMin, upper=rateMax, control=c(fnscale=-1))
					},
		   
		   "medusa" = {
					nodePresent <- phy$edge[,1] %in% phy$edge[,2]
					nodes <- unique(sort(phy$edge[nodePresent,1]))
					nodeDepth <- node.depth(phy)
					nodesByRichness <- cbind(richness= nodeDepth[nodes], node=nodes)
					searchNode <- nodesByRichness[nodesByRichness[,1] >= minCladeSize,2]
					n <- length(y)
					BERateOut <- matrix(NA, nrow=nSplits, ncol=(5+nSplits))
					fullModelOut <- vector(mode="list", length=nSplits)
		   
					for (k in 1:nSplits) {
						cladeMembers <- matrix(NA, ncol=k, nrow=length(phy$edge[,1]))
						if (k==1) {searchNode <- searchNode} else { searchNode <- searchNode[!(searchNode %in% bestNodes)] }
		   
		   
						for (i in searchNode) {
		   
							if (k == 1) { currentNodeIDs <- i } else { currentNodeIDs <- c(bestNodes, i) }
		   
							currentModel <- transformPhylo.ML(y=y, phy=phy, nodeIDs=currentNodeIDs, model="clade") ####
		   
							fullModelOut[[k]] <- rbind(fullModelOut[[k]], c(node=as.integer(i), ML=currentModel$value, currentModel$par))
		   
							currentModel <- list(currentModel, i, cladeMembers)
							currentML <- currentModel[[1]]$value
							AIC <- -2 * currentModel[[1]]$value + 2 * k
							AICc <- -2 * currentModel[[1]]$value + 2 * k + ((2*k*(k+1))/(n-k-1))
		   
		   
							if (i == min(searchNode)) {BERateOut[k, 1:(5+k)] <- c(currentModel[[2]], currentModel[[1]]$value, k+2, AIC, AICc, currentModel[[1]]$par)} else {
		   
							if (currentML > BERateOut[k, 2]) {BERateOut[k, 1:(5+k)] <- c(currentModel[[2]], currentModel[[1]]$value, k+2, AIC, AICc, currentModel[[1]]$par) }
							}
							}
						
						if (k==1) {bestNodes <- BERateOut[k,1]} else {bestNodes <- c(bestNodes, BERateOut[k,1])}
		   
						print(BERateOut)
						}
		   
		   
					BERateSummary <- matrix(NA, ncol=(5+nSplits), nrow=1)
					MLsingle <- likTraitPhylo(y, phy)[[2]]
					AICsingle <- -2 * MLsingle + 2 * k
					AICcsingle <- -2 * MLsingle + 2 * k + ((2*k*(k+1))/(n-k-1)) 
					BERateSummary[1,1:5] <- c(0, MLsingle, 2, AICsingle, AICcsingle)		   		
					BERateSummary <- rbind(BERateSummary, BERateOut)
					colnames(BERateSummary) <- c("node", "ML", "k", "AIC", "AICc",c(BERateSummary[2:nrow(BERateSummary),1]))
		   
					vo <- list(as.data.frame(BERateSummary), fullModelOut)
					}
		   
		   )	
	

	return(vo)
}
			
