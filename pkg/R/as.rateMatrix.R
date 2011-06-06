as.rateMatrix <-
function(phy, x, data) {
	
		# Check that number of taxa matches between data frame and phylogeny
		if(length(phy$tip.label) != length(rownames(data))){stop("Length of phenotypic and of phylogenetic data do not match")}	# check names

	
		# Remove species with NA for the discrete trait from the phylogeny and dataframe
		idxphy <- which(is.na(x) == TRUE)
		if(length(idxphy)>0) { phy <- drop.tip(phy, tip=rownames(data)[idxphy]) } else {phy <- phy }
		
		idxdat <- which(is.na(x) == FALSE)
		data <- data[idxdat,]
	
		# alphabetise tip labels with index return
		speciesPhylo <- sort(phy$tip.label, index.return=TRUE)
	
		# make data frame of discrete trait with species as rownames
		discreteTrait <- data.frame(x, row.names=rownames(data))
	
		# alphabetise species names for trait data with index return
		speciesData <- sort(rownames(discreteTrait), index.return = TRUE) 
	
		if(sum(speciesPhylo$x != speciesData$x) > 0){stop("Taxon names in phylogeny do not match those in data frame")}	# check names

		# data frame of tip states in species alphabetical order
		speState <- data.frame(tip.index = speciesPhylo$ix, data.index = speciesData$ix, row.names=rownames(discreteTrait)[speciesData$ix], tipState=discreteTrait[speciesData$ix, ]) 

		# combines vector of internal nodes with numeric vector of ancestral states taken from node labels excluding root
		ancState <- data.frame(node = c((length(phy$tip.label)+2):max(phy$edge)), ancState = as.numeric(phy$node.label[c(2:length(phy$node.label))]))

		# Put together node and tip indices and ancestral states then reorder
		state<- data.frame(edge.index = c(ancState$node, speState$tip.index), state = c(ancState$ancState, speState$tipState))
		state <- state[sort(state$edge.index, index.return = TRUE)$ix, ]
	
		
		# Put tree edge in node order and then put it all together and order according to the tree edge index
        phyEdge <- sort(phy$edge[,2], index.return = TRUE)
        phyEdgeState <- data.frame(state, phyEdge = phyEdge$x, phyEdgeix = phyEdge$ix)
        phyEdgeState <- phyEdgeState[order(phyEdgeState$phyEdgeix),]
	
		# get numeric vector of branch states (use as factor first that factor levels are scored from 1-i in either numerical of alphabetical order depending on input
		ancestor <- as.numeric(as.factor(phyEdgeState$state)) 
		
		# how many factor levels?
		nlevels <- sort(unique(ancestor))
	
		# set up list to contain vcv matrices for each state
		rateMatrix <- vector(mode="list", length = length(nlevels))
	
		# create state vcv matrices 
		for(i in nlevels) {
				x.anc <- ifelse(ancestor == i, 1, 0)
				state.edge <- x.anc * phy$edge.length
				state.phy <- list(edge=phy$edge, edge.length=state.edge, Nnode=phy$Nnode, tip.label=phy$tip.label)
				class(state.phy) <- "phylo"
			state.matrix <- VCV.array(state.phy)
			class(state.matrix) <- "matrix"
			rateMatrix[[i]] <- state.matrix
		}		
	
	class(rateMatrix) <- "rateMatrix"
	return(rateMatrix)
}

