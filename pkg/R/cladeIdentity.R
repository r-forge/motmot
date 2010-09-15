cladeIdentity <- function(phy, nodeIDs) {
		
		k <- length(nodeIDs)
		cladeMembers <- matrix(NA, ncol=k, nrow=length(phy$edge[,1]))

		for (i in 1:k) {
			nodeShiftID <- c(nodeIDs[i], node.descendents(x=nodeIDs[i], phy=phy))
			cladeMembers[,i] <- as.numeric(phy$edge[,2] %in% nodeShiftID)
						}
    	
    	originalOrder <- apply(cladeMembers, 2, sum)
    	richnessOrder <- sort(originalOrder, decreasing=FALSE, index.return=TRUE)

    	cladeMembersOrdered <- matrix(cladeMembers[, richnessOrder$ix], ncol=length(nodeIDs))
    
    	for (i in 1:k) {		
			for (j in 1:k) {
				if (i!=j) {nestedClade <- which(cladeMembersOrdered[,i]&cladeMembersOrdered[, j]==1)    
					cladeMembersOrdered[nestedClade, j] <- 0
				}}}	
				
			 cladeMembers <- cladeMembersOrdered[, sort(richnessOrder$ix, index.return=TRUE)$ix]
			 cladeMembers <- 	matrix(cladeMembers, ncol=length(nodeIDs))		
			return(cladeMembers)
			}
			
