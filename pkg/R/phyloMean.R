phyloMean <-
function(rateData, rate=NULL, common.mean=FALSE) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
				
		if(length(rate) != length(rateData$Vmat)){stop("The number of rates defined differs from the number of rate matrices")}
	
		y <- rateData$y
		x <- as.factor(rateData$x)

		V <- transformRateMatrix(rateData, rate)
		x <- make.anc(y, x)
			
		if(common.mean==FALSE) {x <- x} else { x <- rep(1, length(x[,1]))}

			iV <- solve(V)
			xVix <- crossprod(x, iV %*% x)
			xViy <- crossprod(x, iV %*% y)
			mu <- solve(xVix) %*% xViy 
			return(mu)
		}

