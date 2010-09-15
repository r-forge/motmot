transformPhylo.ML <- function(y, phy, model=NULL, modelCIs=TRUE, nodeIDs=NULL, rateType=NULL, minCladeSize=1, nSplits=10, restrictNode=NULL, lowerBound=NULL, upperBound=NULL){
	
	
	bounds <- matrix(c(1e-7,1,1e-7,1,1e-7,5,1e-7,10,0,1,1e-7,200), 6, 2, byrow=TRUE)
	rownames(bounds) <- c("kappa", "lambda", "delta", "alpha", "psi", "rate")
	
	switch(model,
		   
		   "bm" = {
					out <- likTraitPhylo(y, phy)
					names(out) <- c("brownianVariance","logLikelihood")
					},
		   
		   "kappa" = {
					kappa <- 1
					if (is.null(lowerBound)) { lowerBound <- bounds["kappa", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["kappa", 2] }
					var.funkappa <- function(kappa) { return(transformPhylo.ll(y, phy, kappa, model="kappa")[[2]])}
					vo <- optimise(var.funkappa, interval=c(lowerBound, upperBound), maximum=TRUE)
		   
					if (modelCIs==TRUE) {
		   
						foo <- function(param) {
							ll <- transformPhylo.ll(y, phy, model="kappa", kappa=param)$logLikelihood
							return(ll - vo$objective + 1.92)
							}
		   
						if(foo(lowerBound) < 0) { 
							LCI <- uniroot(foo, interval = c(lowerBound, vo$maximum))$root 
							} else { LCI <- NA }
						if(foo(upperBound) < 0) {
							UCI <- uniroot(foo, interval = c(vo$maximum, upperBound))$root
						} else { UCI <- NA }
							   
					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Kappa")
		   
					out$Kappa <- matrix(NA, 1, 3, byrow=TRUE)
					colnames(out$Kappa) <- c("MLKappa", "LowerCI", "UpperCI")
		   
					out$MaximumLikelihood <- vo$objective
					out$Kappa[1,] <- c(vo$maximum, LCI, UCI)
		   
		   
					if (any(is.na(out$Kappa[1,2:3]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}

					} else { 

					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Kappa")
		   
					out$Kappa <- matrix(NA, 1, 1, byrow=TRUE)
					colnames(out$Kappa) <- c("MLKappa")
		   
					out$MaximumLikelihood <- vo$objective
					out$Kappa[1,] <- c(vo$maximum)
		   
					}
					},

		   "lambda" = {
					lambda <- 1
					if (is.null(lowerBound)) { lowerBound <- bounds["lambda", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["lambda", 2] }
					var.funlambda <- function(lambda) { return(transformPhylo.ll(y=y, phy=phy, lambda=lambda, model="lambda")[[2]])}
					vo <- optimise(var.funlambda, interval=bounds["lambda",], maximum=TRUE)
		   
					if (modelCIs==TRUE) {
		   
						foo <- function(param) {
							ll <- transformPhylo.ll(y, phy, model="lambda", lambda=param)$logLikelihood
							return(ll - vo$objective + 1.92)
							}
		   
						if(foo(lowerBound) < 0) { 
							LCI <- uniroot(foo, interval = c(lowerBound, vo$maximum))$root 
							} else { LCI <- NA }
						if(foo(upperBound) < 0) {
							UCI <- uniroot(foo, interval = c(vo$maximum, upperBound))$root
							} else { UCI <- NA }
		   
					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Lambda")
		   
					out$Lambda <- matrix(NA, 1, 3, byrow=TRUE)
					colnames(out$Lambda) <- c("MLLambda", "LowerCI", "UpperCI")
		   
					out$MaximumLikelihood <- vo$objective
					out$Lambda[1,] <- c(vo$maximum, LCI, UCI)
		   
					
					if (any(is.na(out$Lambda[1,2:3]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}

					} else { 

					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Lambda")
		   
					out$Lambda <- matrix(NA, 1, 1, byrow=TRUE)
					colnames(out$Lambda) <- c("MLLambda")
		   
					out$MaximumLikelihood <- vo$objective
					out$Lambda[1,] <- c(vo$maximum)
		   
					}
					},
		   
		   "delta" = {
					delta <- 1
					if (is.null(lowerBound)) { lowerBound <- bounds["delta", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["delta", 2] }
					var.fundelta <- function(delta) { return(transformPhylo.ll(y=y, phy=phy, delta=delta, model="delta")[[2]])}
					vo <- optimise(var.fundelta, interval=bounds["delta",], maximum=TRUE)
		   
					if (modelCIs==TRUE) {
		   
						foo <- function(param) {
							ll <- transformPhylo.ll(y, phy, model="delta", delta=param)$logLikelihood
							return(ll - vo$objective + 1.92)
							}
		   
						if(foo(lowerBound) < 0) { 
							LCI <- uniroot(foo, interval = c(lowerBound, vo$maximum))$root 
							} else { LCI <- NA }
						if(foo(upperBound) < 0) {
							UCI <- uniroot(foo, interval = c(vo$maximum, upperBound))$root
							} else { UCI <- NA }
		   
					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Delta")
		   
					out$Delta <- matrix(NA, 1, 3, byrow=TRUE)
					colnames(out$Delta) <- c("MLDelta", "LowerCI", "UpperCI")
		   
					out$MaximumLikelihood <- vo$objective
					out$Delta[1,] <- c(vo$maximum, LCI, UCI)
		   
		   
					if (any(is.na(out$Delta[1,2:3]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}

					} else { 

					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Delta")
		   
					out$Delta <- matrix(NA, 1, 1, byrow=TRUE)
					colnames(out$Delta) <- c("MLDelta")
		   
					out$MaximumLikelihood <- vo$objective
					out$Delta[1,] <- c(vo$maximum)
		   
					}
					},
		   
		   "OU" = {
					alpha <- 0.01
					if (is.null(lowerBound)) { lowerBound <- bounds["alpha", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["alpha", 2] }
					var.funOU <- function(alpha) { return(transformPhylo.ll(y=y, phy=phy, alpha=alpha, model="OU")[[2]])}
					vo <- optimise(var.funOU, interval=bounds["alpha",], maximum=TRUE)
		   
					if (modelCIs==TRUE) {
		   
						foo <- function(param) {
							ll <- transformPhylo.ll(y, phy, model="OU", alpha=param)$logLikelihood
							return(ll - vo$objective + 1.92)
							}
		   
						if(foo(lowerBound) < 0) { 
							LCI <- uniroot(foo, interval = c(lowerBound, vo$maximum))$root 
							} else { LCI <- NA }
						if(foo(upperBound) < 0) {
							UCI <- uniroot(foo, interval = c(vo$maximum, upperBound))$root
							} else { UCI <- NA }

		   
					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Alpha")
		   
					out$Alpha <- matrix(NA, 1, 3, byrow=TRUE)
					colnames(out$Alpha) <- c("MLAlpha", "LowerCI", "UpperCI")
		   
					out$MaximumLikelihood <- vo$objective
					out$Alpha[1,] <- c(vo$maximum, LCI, UCI)
		   
		   
					if (any(is.na(out$Alpha[1,2:3]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}

					} else { 

					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Alpha")
		   
					out$Alpha <- matrix(NA, 1, 1, byrow=TRUE)
					colnames(out$Alpha) <- c("MLAlpha")
		   
					out$MaximumLikelihood <- vo$objective
					out$Alpha[1,] <- c(vo$maximum)
		   
					}
					},
		   
		   "psi" = {
					psi <- 1
					if (is.null(lowerBound)) { lowerBound <- bounds["psi", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["psi", 2] }
					var.funpsi <- function(psi) { return(transformPhylo.ll(y=y, phy=phy, psi=psi, model="psi")[[2]])}
					vo <- optimise(var.funpsi, interval=bounds["psi",], maximum=TRUE)
		   
					if (modelCIs==TRUE) {
		   
					foo <- function(param) {
						ll <- transformPhylo.ll(y, phy, model="psi", psi=param)$logLikelihood
						return(ll - vo$objective + 1.92)
						}
		   
				if(foo(lowerBound) < 0) { 
					LCI <- uniroot(foo, interval = c(lowerBound, vo$maximum))$root 
					} else { LCI <- NA }
				if(foo(upperBound) < 0) {
					UCI <- uniroot(foo, interval = c(vo$maximum, upperBound))$root
					} else { UCI <- NA }
		   
				out <- vector(mode="list", length=2)
				names(out) <- c("MaximumLikelihood", "psi")
		   
				out$psi <- matrix(NA, 1, 3, byrow=TRUE)
				colnames(out$psi) <- c("MLpsi", "LowerCI", "UpperCI")
		   
				out$MaximumLikelihood <- vo$objective
				out$psi[1,] <- c(vo$maximum, LCI, UCI)
		   
		   
				if (any(is.na(out$psi[1,2:3]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}
		   
				} else { 
		   
				out <- vector(mode="list", length=2)
				names(out) <- c("MaximumLikelihood", "psi")
		   
				out$psi <- matrix(NA, 1, 1, byrow=TRUE)
				colnames(out$psi) <- c("MLpsi")
		   
				out$MaximumLikelihood <- vo$objective
				out$psi[1,] <- c(vo$maximum)
		   
				}
				},
		   
		   "free" = {
					if (is.null(lowerBound)) { lowerBound <- bounds["rate", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["rate", 2] }
					branchRates <- rep(1,length(phy$edge.length))
					var.funfree <- function(branchRates) { return(transformPhylo.ll(y, phy, branchRates=branchRates, model="free")[[2]])}
					out <- optim(branchRates, var.funfree, method="L-BFGS-B", lower=lowerBound, upper=upperBound, control=c(fnscale=-1))
					},
		   
		   "clade" = {
					if (is.null(lowerBound)) { lowerBound <- bounds["rate", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["rate", 2] }
					cladeRates <- rep(1,length(nodeIDs))
					var.funclade <- function(cladeRates) { return(transformPhylo.ll(y, phy, nodeIDs=nodeIDs, cladeRates=cladeRates, model="clade")[[2]])}
					vo <- optim(cladeRates, var.funclade, method="L-BFGS-B", lower=lowerBound, upper=upperBound, control=c(fnscale=-1))
					
		   		   
					if (modelCIs==TRUE) {
		   
					foo <- function(param) {
						ll <- transformPhylo.ll(y, phyClade, model="clade", nodeIDs=SingleNode, rateType=whichRateType, cladeRates=param)$logLikelihood
						return(ll - vo$value + 1.92)
						}
		   
					phyClade <- transformPhylo(phy, model="clade", nodeIDs=nodeIDs, rateType=rateType, cladeRates=vo$par)
		   
					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Rates")
		   
					out$Rates <- matrix(NA, length(nodeIDs), 4, byrow=TRUE)
					colnames(out$Rates) <- c("node", "MLRate", "LowerCI", "UpperCI")
		   
					out$MaximumLikelihood <- vo$value
		   
						for (i in 1:length(nodeIDs)) {
		   
							SingleNode <- nodeIDs[i]
							whichRateType <- rateType[i]
							phyClade <- transformPhylo(phyClade, model="clade", nodeIDs=SingleNode, rateType=whichRateType, cladeRates=1/vo$par[i])
				
							if(foo(lowerBound) < 0) { 
								LCI <- uniroot(foo, interval = c(lowerBound, vo$par[i]))$root 
								} else { LCI <- NA }
							if(foo(upperBound) < 0) {
								UCI <- uniroot(foo, interval = c(vo$par[i], upperBound))$root
								} else { UCI <- NA }
					
							out$Rates[i,] <- c(nodeIDs[i], vo$par[i], LCI, UCI)
							}
					
					if (any(is.na(out$Rates[,3:4]))) { warning("Confidence limits fall outside the current parameter bounds - consider changing lowerBound and/or upperBound")}
		   
					} else { 
					out <- vector(mode="list", length=2)
					names(out) <- c("MaximumLikelihood", "Rates")
		   
					out$Rates <- matrix(NA, length(nodeIDs), 2, byrow=TRUE)
					colnames(out$Rates) <- c("node", "MLRate")
		   
					out$MaximumLikelihood <- vo$value
					out$Rates[,1] <- nodeIDs
					out$Rates[,2] <- vo$par
					}
					},
		   
		   "tm2" = {
					if (is.null(lowerBound)) { lowerBound <- bounds["rate", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["rate", 2] }
			
					nodes <- unique(sort(phy$edge[,2]))
					nodeDepth <- node.depth(phy)
					nodesByRichness <- cbind(richness= nodeDepth[nodes], node=nodes)
					searchNode <- nodesByRichness[nodesByRichness[,1] >= minCladeSize,2]
		   
		   
						if (is.null(restrictNode)==FALSE) {
		   
							ngroups <- length(restrictNode)
							cm <- clade.matrix(phy)
		   
							colnames(cm$clade.matrix) <- cm$tip.label
							cmMat <- cm$clade.matrix
		   
							skipNodes <- numeric()
		   
							for (i in 1:ngroups) {
								matDrop <- cmMat[,restrictNode[[i]]] 
								nodeDrop <- rowSums(matDrop) < length(restrictNode[[i]]) & rowSums(matDrop) > 1
								skipNodes <- c(skipNodes, as.numeric(rownames(matDrop)[nodeDrop]))
								}
						
							searchNode <- setdiff(searchNode, skipNodes)
							}
		   
		   
					n <- length(y)
					BERateOut <- matrix(NA, nrow=nSplits, ncol=(6+nSplits))
					fullModelOut <- vector(mode="list", length=nSplits)
		   
					for (k in 1:nSplits) {
						cladeMembers <- matrix(NA, ncol=k, nrow=length(phy$edge[,1]))
						if (k==1) {searchNode <- searchNode} else { searchNode <- searchNode[!(searchNode %in% bestNodes)] }
		   
		   
						for (i in searchNode) {
		   
							if (k == 1) { currentNodeIDs <- i } else { currentNodeIDs <- c(bestNodes, i) }
		   
							cladeRates <- rep(1,length(currentNodeIDs))
							var.funclade <- function(cladeRates) { return(transformPhylo.ll(y, phy, nodeIDs=currentNodeIDs, rateType=rep("clade", length(currentNodeIDs)), cladeRates=cladeRates, model="clade")[[2]])}
							var.funbranch <- function(cladeRates) { return(transformPhylo.ll(y, phy, nodeIDs=currentNodeIDs, rateType=rep("branch", length(currentNodeIDs)), cladeRates=cladeRates, model="clade")[[2]])}

							currentCladeModel <- optim(cladeRates, var.funclade, method="L-BFGS-B", lower=lowerBound, upper=upperBound, control=c(fnscale=-1))
							currentBranchModel <- optim(cladeRates, var.funbranch, method="L-BFGS-B", lower=lowerBound, upper=upperBound, control=c(fnscale=-1))

		   
							fullModelOut[[k]] <- rbind(fullModelOut[[k]], c(node=as.integer(i), shiftPos=1, ML= currentCladeModel$value, currentCladeModel$par))
							fullModelOut[[k]] <- rbind(fullModelOut[[k]], c(node=as.integer(i), shiftPos=2, ML= currentBranchModel$value, currentBranchModel$par))
		   
							if (currentCladeModel$value > currentBranchModel$value) { currentModel <- currentCladeModel; shiftPos=1 }
							if (currentCladeModel$value < currentBranchModel$value) { currentModel <- currentBranchModel; shiftPos=2 }

							param <- k+2
		   
							currentModel <- list(currentModel, i, cladeMembers)
							currentML <- currentModel[[1]]$value
							AIC <- -2 * currentModel[[1]]$value + 2 * param
							AICc <- -2 * currentModel[[1]]$value + 2 * param + ((2*param*(param+1))/(n-param-1))
		   
		   
							if (i == min(searchNode)) {BERateOut[k, 1:(6+k)] <- c(currentModel[[2]], shiftPos, currentModel[[1]]$value, k+2, AIC, AICc, currentModel[[1]]$par)} else {
		   
							if (currentML > BERateOut[k, 3]) {BERateOut[k, 1:(6+k)] <- c(currentModel[[2]], shiftPos, currentModel[[1]]$value, k+2, AIC, AICc, currentModel[[1]]$par) }
							}
							}
						
						if (k==1) {bestNodes <- BERateOut[k,1]} else {bestNodes <- c(bestNodes, BERateOut[k,1])}
		   
						print(BERateOut)
						}
		   
		   
					BERateSummary <- matrix(NA, ncol=(6+nSplits), nrow=1)
					BERateSummary <- as.data.frame(BERateSummary)
					MLsingle <- likTraitPhylo(y, phy)[[2]]
					AICsingle <- -2 * MLsingle + 2 * 2
					AICcsingle <- -2 * MLsingle + 2 * 2 + ((2*2*(2+1))/(n-2-1)) 
					BERateSummary[1,1:6] <- c(0, 1, MLsingle, 2, AICsingle, AICcsingle)		   		
					BERateSummary <- rbind(BERateSummary, BERateOut)
					colnames(BERateSummary) <- c("node", "shiftPos", "ML", "k", "AIC", "AICc",c(BERateSummary[2:nrow(BERateSummary),1]))
		   
		   
					if (sum(BERateSummary[,"shiftPos"]==1) > 0) {BERateSummary[which(BERateSummary[,"shiftPos"]==1), "shiftPos"] <- "clade"}
					if (sum(BERateSummary[,"shiftPos"]==2) > 0) {BERateSummary[which(BERateSummary[,"shiftPos"]==2), "shiftPos"] <- "branch"}

					out <- list(as.data.frame(BERateSummary), fullModelOut, y, phy)
					},

		   
		   "tm1" = {
					if (is.null(lowerBound)) { lowerBound <- bounds["rate", 1] }
					if (is.null(upperBound)) { upperBound <- bounds["rate", 2] }
		   
					nodes <- unique(sort(phy$edge[,2]))
					nodeDepth <- node.depth(phy)
					nodesByRichness <- cbind(richness= nodeDepth[nodes], node=nodes)
					searchNode <- nodesByRichness[nodesByRichness[,1] >= minCladeSize,2]
		   
		   
					if (is.null(restrictNode)==FALSE) {
		   
						ngroups <- length(restrictNode)
						cm <- clade.matrix(phy)
		   
						colnames(cm$clade.matrix) <- cm$tip.label
						cmMat <- cm$clade.matrix
		   
						skipNodes <- numeric()
			
						for (i in 1:ngroups) {
							matDrop <- cmMat[,restrictNode[[i]]] 
							nodeDrop <- rowSums(matDrop) < length(restrictNode[[i]]) & rowSums(matDrop) > 1
							skipNodes <- c(skipNodes, as.numeric(rownames(matDrop)[nodeDrop]))
							}
		   
						searchNode <- setdiff(searchNode, skipNodes)
						}
		   
		   
					n <- length(y)
					BERateOut <- matrix(NA, nrow=nSplits, ncol=(6+nSplits))
					fullModelOut <- vector(mode="list", length=nSplits)
		   
					for (k in 1:nSplits) {
						cladeMembers <- matrix(NA, ncol=k, nrow=length(phy$edge[,1]))
						if (k==1) {searchNode <- searchNode} else { searchNode <- searchNode[!(searchNode %in% bestNodes)] }
		   
		   
						for (i in searchNode) {
		   
							if (k == 1) { currentNodeIDs <- i } else { currentNodeIDs <- c(bestNodes, i) }
		   
							cladeRates <- rep(1,length(currentNodeIDs))
							var.funclade <- function(cladeRates) { return(transformPhylo.ll(y, phy, nodeIDs=currentNodeIDs, rateType=rep("clade", length(currentNodeIDs)), cladeRates=cladeRates, model="clade")[[2]])}
		   
							currentCladeModel <- optim(cladeRates, var.funclade, method="L-BFGS-B", lower=lowerBound, upper=upperBound, control=c(fnscale=-1))
		   
		   
							fullModelOut[[k]] <- rbind(fullModelOut[[k]], c(node=as.integer(i), shiftPos=1, ML= currentCladeModel$value, currentCladeModel$par))
		   
							currentModel <- currentCladeModel 
							shiftPos=1
		   
							param <- k+2
		   
							currentModel <- list(currentModel, i, cladeMembers)
							currentML <- currentModel[[1]]$value
							AIC <- -2 * currentModel[[1]]$value + 2 * param
							AICc <- -2 * currentModel[[1]]$value + 2 * param + ((2*param*(param+1))/(n-param-1))
		   
		   
							if (i == min(searchNode)) {BERateOut[k, 1:(6+k)] <- c(currentModel[[2]], shiftPos, currentModel[[1]]$value, k+2, AIC, AICc, currentModel[[1]]$par)} else {
		   
							if (currentML > BERateOut[k, 3]) {BERateOut[k, 1:(6+k)] <- c(currentModel[[2]], shiftPos, currentModel[[1]]$value, k+2, AIC, AICc, currentModel[[1]]$par) }
						}
						}
		   
		   
					if (k==1) {bestNodes <- BERateOut[k,1]} else {bestNodes <- c(bestNodes, BERateOut[k,1])}
		   
					print(BERateOut)
					}
		   
		   
				BERateSummary <- matrix(NA, ncol=(6+nSplits), nrow=1)
				BERateSummary <- as.data.frame(BERateSummary)
				MLsingle <- likTraitPhylo(y, phy)[[2]]
				AICsingle <- -2 * MLsingle + 2 * 2
				AICcsingle <- -2 * MLsingle + 2 * 2 + ((2*2*(2+1))/(n-2-1)) 
				BERateSummary[1,1:6] <- c(0, 1, MLsingle, 2, AICsingle, AICcsingle)		   		
				BERateSummary <- rbind(BERateSummary, BERateOut)
				colnames(BERateSummary) <- c("node", "shiftPos", "ML", "k", "AIC", "AICc",c(BERateSummary[2:nrow(BERateSummary),1]))
		   
		   
				if (sum(BERateSummary[,"shiftPos"]==1) > 0) {BERateSummary[which(BERateSummary[,"shiftPos"]==1), "shiftPos"] <- "clade"}
		   
				out <- list(as.data.frame(BERateSummary), fullModelOut, y, phy)
				}
		   
		   
			)	
	
	return(out)
}








