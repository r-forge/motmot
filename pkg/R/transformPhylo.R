transformPhylo <- function(phy, model=NULL, kappa=NULL, lambda=NULL, delta=NULL, alpha=NULL, psi=NULL, nodeIDs=NULL, rateType=NULL, branchRates=NULL, cladeRates=NULL){
	
	switch(model,
		   "kappa" = {
					 phy$edge.length <- phy$edge.length ^ kappa
					 },
		   		   
		   "lambda" = {
					
					if (is.ultrametric(phy)) {
						rootOrig <- max(branching.times(phy))
						tips <- match(c(1:Ntip(phy)), phy$edge[,2])
						phy$edge.length <- phy$edge.length * lambda
						phy$edge.length[tips] <- phy$edge.length[tips] + (rootOrig * (1-lambda))
						}  
		   
					if (is.ultrametric(phy)==FALSE)	{
						tips <- match(c(1:Ntip(phy)), phy$edge[,2])
						cladeMat <- clade.matrix(phy)
						branchHeights <- rep(NA, Ntip(phy))
						for (i in 1:Ntip(phy)) { 	branchHeights[i] <- sum(cladeMat$edge.length[cladeMat$clade.matrix[,i]==1]) }
		   
						phy$edge.length <- phy$edge.length * lambda
						phy$edge.length[tips] <- phy$edge.length[tips] + (branchHeights * (1-lambda))
						}
					},
		   
		   "delta" = {
					times <- branching.times(phy)
					times <- max(times) - times
					tips <- length(phy$tip.label)
					res <- phy
		 
					for (i in 1:length(phy$edge.length)) {
							bl <- phy$edge.length[i]
							age <- times[phy$edge[i, 1] - tips]
							res$edge.length[i] <- (age + bl)^delta - age^delta
							}
					phy <- res
					},
					
		   
		   "free" = {
					branchRates <- branchRates + (1-min(branchRates))
					phy$edge.length <- phy$edge.length * branchRates
					},
		     
		   "clade" = {
					cladeMembers <- cladeIdentity(phy=phy, nodeIDs=nodeIDs)
					if (is.null(rateType)) { rateType <- rep("clade", length(nodeIDs))} else {rateType <- rateType}
		   
					for (i in 1:length(cladeRates)) {
						if (rateType[i]=="clade") {phy$edge.length[cladeMembers[,i]==1] <- phy$edge.length[cladeMembers[,i]==1] * cladeRates[i]}
						if (rateType[i]=="branch"){phy$edge.length[phy$edge[,2]==nodeIDs[i]] <- phy$edge.length[phy$edge[,2]==nodeIDs[i]] * cladeRates[i]}
						}
					},
		   
		   
		   "OU" = {
					times <- branching.times(phy)
					names(times) <- (Ntip(phy) + 1):(Ntip(phy) + Nnode(phy))
					Tmax<-times[1]
					phy2<-phy
		   
					for (i in 1:length(phy$edge.length)) {
						bl <- phy$edge.length[i]
						age <- times[which(names(times) == phy$edge[i, 1])]
						t1 <- max(times) - age
						t2 <- t1+bl
					phy2$edge.length[i] <- (1/(2*alpha))*exp(-2*alpha * (Tmax-t2)) * (1 - exp(-2 * alpha * t2)) - 
					(1/(2*alpha))*exp(-2*alpha * (Tmax-t1)) * (1 - exp(-2 * alpha * t1))
					}
					phy <- phy2
					},
		   
		   "psi" = {
					bt <- branching.times(phy)
					bdrates <- bd(bt)
					mu <- bdrates$r1 * bdrates$a	# estimate of extinction rate from the tree
					lambda <- bdrates$r1 - mu #  estimate of speciation rate from the tree
					probNoDescendents <- mu*(exp((lambda-mu)*phy$edge.length)-1) / ((lambda*exp((lambda-mu)*phy$edge.length))-mu) # probability that a branch has no descendents
					sh <- (lambda*probNoDescendents) * phy$edge.length
					phy2 <- phy
					phy2$edge.length <- psi * (phy$edge.length^0 + sh) + (1-psi)*phy$edge.length
					phy <- phy2
					}		   
		   
	   
		   )
		   
		   return(phy)
		   
		   }
			
