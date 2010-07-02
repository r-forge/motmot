fairProportions <- function (phy) {
	
	treeMatrix <- clade.matrix(phy)
	fpEdgeVals <- treeMatrix$edge.length / apply(treeMatrix$clade.matrix, 1, sum) 
	fpEdgeMatrix <- fpEdgeVals * treeMatrix$clade.matrix
	fpTips <- apply(fpEdgeMatrix, 2, sum)
	
	names(fpTips) <- phy$tip.label
	return(fpTips)
	
}