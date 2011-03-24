make.likRatePhylo <-
function(rateData, fixed, common.mean=FALSE) {
								
		op <- fixed 

		function(rate){
			op[!fixed] <- rate[!fixed]
			
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
			
			V <- rateData$Vmat
			v1 <- V[[1]]
			nV <- length(op)
			
			rateMats <- vector(mode="list", length = nV)
			retMat <- matrix(0, nrow = dim(v1)[1], ncol = dim(v1)[2])
			
			for(i in 1:nV) {
				rateMats[[i]] <- op[i] * V[[i]]  
				retMat <- retMat + rateMats[[i]]
			}
			
			meserr <- rateData$meserr
			if(is.null(meserr) == FALSE) { 
				y <- rateData$y
				diag(retMat) <- diag(retMat) + (meserr^2) / (var(y)/max(retMat)) }
			
			x <- make.anc(y, x)
			
			if(common.mean==FALSE) {
				x <- x} else { x <- rep(1, length(x[,1]))}
						
			logDetV <- logDetFun(retMat)
			iV <- solve(retMat)
			xVix <- crossprod(x, iV %*% x)
			xViy <- crossprod(x, iV %*% y)
			mu <- solve(xVix) %*% xViy 

			e <- y - x %*% mu
			s2 <- crossprod(e, iV %*% e)
			n <- length(y) 
			phylo.var <- ( s2 / (n - k) )
		
			n <- length(y)
			ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
			
			ypred <- x%*%mu
			lik.RatePhylo <- ( list(ll = ll, mu = mu, phylo.var = phylo.var) )
			return(-1 * lik.RatePhylo$ll)			
	}}

