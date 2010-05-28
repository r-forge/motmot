dummy.RatePhylo <-
function(rateData, rate=NULL, group.means=NULL) {

	if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		
	V <- transformRateMatrix(rateData, rate=rate)
	expect.sd <- sqrt(mean(V[upper.tri(V)]))
			
	if (is.null(group.means)) {ydum <- as.matrix(t(rmvnorm(1, sigma =  (V) ))) } 
	
		else {
			x.means <- unique(rateData$x)
			n.means <- length(x.means)
			samp.means <- rep(NA, length(rateData$x))
			ydum <- vector(mode="list", length=length(group.means))
			for (i in 1:n.means) {
				samp.means[which(rateData$x == (i-1))] <- rep(0+(expect.sd*group.means[i]), length(which(rateData$x == (i-1))))
				}	
			
		ydum <- as.matrix(t(rmvnorm(1, mean=samp.means, sigma =  (V) ))) 		}	
		return(ydum)
	}

