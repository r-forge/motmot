transformPhylo <- function(phy, model=NULL, kappa=NULL, lambda=NULL, delta=NULL, nodeIDs=NULL, branchRates=NULL, cladeRates=NULL){
	
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
					for (i in 1:length(cladeRates)) {
					phy$edge.length[cladeMembers[,i]==1] <- phy$edge.length[cladeMembers[,i]==1] * cladeRates[i]}
					}
		   )
		   
		   return(phy)
		   
		   }
			
