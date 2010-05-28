phyloVar <-
function(rateData, rate=NULL, common.mean=FALSE) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(length(rate) != length(rateData$Vmat)){stop("The number of rates defined differs from the number of rate matrices")}
			
		y <- rateData$y
		x <- as.factor(rateData$x)
						
		if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
			
		V <- transformRateMatrix(rateData, rate)
		x <- make.anc(y, x)
			
		if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}
		
		mu <- phyloMean(rateData, rate, common.mean=common.mean)
			
			iV <- solve(V)
			e <- y - x %*% mu
			s2 <- crossprod(e, iV %*% e)
			n <- length(y) 
			phylo.var <- ( s2 / (n - k) )
			return(phylo.var)
			}

