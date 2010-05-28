optim.likRatePhylo <-
function(rateData, rate=NULL, fixed = NULL, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE) {
		
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(is.null(fixed))  { op <- fixed <- c(rep(FALSE,length(rateData$Vmat) - 1), FALSE) } else { op <- fixed }	
			
		mvl <- make.likRatePhylo(rateData, op, common.mean=common.mean)

			vo <- try(optim(rate, mvl, method = "L-BFGS-B", lower = rateMIN, upper = rateMAX))
			MLRate <- vo$par
			
			fixed[which(fixed==FALSE)] <- MLRate[which(fixed==FALSE)]
			MLRate <- fixed

			ML <- -vo$value
			convergence <- vo$convergence
			n <- length(rateData$y)
			
				if(length(op)!=length(which(op==FALSE))) {
					if(common.mean==TRUE) {k <- 2 + length(which(op==FALSE))
							} else { k <- (length(which(op==FALSE)) +1 + length(op)) }
							
							} else {
					
					if(common.mean==TRUE) {k <- 1 + length(which(op==FALSE))
							} else { k <- (length(which(op==FALSE)) + length(op)) }
							}
							
			aic <- -2 * ML + 2 * k
			aicc <- -2 * ML + 2 * k + ((2*k*(k+1))/(n-k-1))
			
			
	ML.RatePhylo <- list(MLRate = MLRate, Max.lik = ML, aic = aic, aicc = aicc, convergence=convergence, n.parameters = k)
	return(ML.RatePhylo)
	
}

