traitMedusaSummary <- function (phy=NULL, traitMedusaObject=NULL, cutoff=4, AICc=TRUE, plotTree=TRUE, stretchTree=FALSE) {
		
	breaks <- numeric()

    if (AICc==TRUE) { datCol = 5 } else { datCol = 4 }
    
	for (i in 2:dim(traitMedusaObject[[1]])[1]) {
        if ((traitMedusaObject[[1]][i-1, datCol] - traitMedusaObject[[1]][i, datCol]) < cutoff) 
            break
        bestModel <- traitMedusaObject[[1]][i,]
        
    }
        
    nSplit <- as.integer(rownames(bestModel))-1
    optimalModelID <- which(traitMedusaObject[[2]][[nSplit]][,2]==max(traitMedusaObject[[2]][[nSplit]][,2]))
    optimalModel <- traitMedusaObject[[2]][[nSplit]][optimalModelID,]
    optimalModel <- matrix(optimalModel, nrow=1)
       
    optimalShiftNodes <- traitMedusaObject[[1]][2:(nSplit+1),1]
    

     
    cladeRates <- optimalModel[,3:(2+nSplit)]
    bestModelOut <- bestModel[,2:(5+ nSplit)]
    rownames(bestModelOut) <- "BestModel"
    optimalTree <- transformPhylo(phy=phy, nodeIDs=optimalShiftNodes, cladeRates=cladeRates, model="clade")
    
    
    
    rateColours <- rev(rainbow(1000)[1:700])
    cladeMembers <- cladeIdentity(phy, optimalShiftNodes)
    edgeColours <- rep("black", length(phy$edge[,1]))
    
    RelRateLimit <- max(c(max(cladeRates), 1/min(cladeRates)))
   		
   		    
    for (i in 1:length(cladeRates)) {
    		
    		if (length(cladeRates) == 1) { if (cladeRates[i] < 1) {colourID <- 1} else {colourID <- 700} }
    			
    			 else { colourID <- round(350 + log10(cladeRates[i]) * (350/log10(RelRateLimit)))} 
    			 
    			if (colourID==0) { colourID <- 1 }
    			if (colourID==701) { colourID <- 700 }
    	
			edgeColours[cladeMembers[,i]==1] <- rateColours[colourID] 
						}
    
    
    if (plotTree) { 
    	if (stretchTree) { plot(optimalTree, edge.color=edgeColours, edge.width=2) } else { plot(phy, edge.color=edgeColours, edge.width=2) } 
    	}
		
		return(list(bestModelOut, optimalTree))
		
		 }
