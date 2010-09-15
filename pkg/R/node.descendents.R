node.descendents <- function (x, phy, tip.labels=FALSE) 
{
    NallNodes <- max(phy$edge)
    Ntips <- max(phy$edge) - phy$Nnode
    nodes <- numeric()
    while (length(intN <- x[x > Ntips]) > 0) {
        minIntN <- min(intN)
        childOfMinIntN <- with(phy, edge[, 2][which(edge[, 1] == minIntN)])
        x <- c(x[x != minIntN], childOfMinIntN)
        nodes <- c(nodes, childOfMinIntN)
    }
	
    nodes <- unique(nodes)
	
	if (tip.labels==TRUE) {	
		
		nodesTips <- vector(mode="list", length=2)
		nodesTips[[1]] <- nodes[nodes > Ntips]
		nodesTips[[2]] <- with(phy, tip.label[nodes[nodes <= Ntips]])
		return(nodesTips)		
		
	} else { return(nodes) }
	
    
}
