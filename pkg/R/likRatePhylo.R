likRatePhylo <-
function(rateData, rate=NULL, common.mean=FALSE) {
	
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	

			
			logDetFun <- function(mat) {
			svdMat <- La.svd(mat)
			d <- svdMat$d
			n <- length(d)
			logDet <- sum(log(d))
			return(logDet)
			}
			y <- rateData$y
			x <- as.factor(rateData$x)
			if (common.mean==FALSE) {k <- nlevels(x)} else {k <- 1}
			V <- transformRateMatrix(rateData, rate)
			x <- make.anc(y, x)
						
			logDetV <- logDetFun(V)
			mu <- phyloMean(rateData, rate, common.mean)
			s2 <- phyloVar(rateData, rate, common.mean)
			n <- length(x[,1])
			ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(s2) - logDetV / 2.0 - (n - k)/2.0
			
			lik.RatePhylo <- ( list(ll = ll, mu = mu, s2 = s2) )
			return(lik.RatePhylo)
	}

