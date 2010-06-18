node.descendents <- function (x, phy) 
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
    return(nodes)
    }