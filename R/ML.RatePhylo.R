ML.RatePhylo <-
function(rateData, rate=NULL, fixed = NULL, pretty = TRUE, rateMIN = 0.001, rateMAX = 50, common.mean=FALSE) {
						
		if(is.null(rate))  { rate <- c(rep(1,length(rateData$Vmat))) } else { rate <- rate }	
	
		if(is.null(fixed))  { op <- fixed <- c(rep(FALSE,length(rateData$Vmat) - 1), TRUE) } else { op <- fixed }		
		
				ovp <- optim.likRatePhylo(rateData, rate, fixed, rateMIN, rateMAX, common.mean=common.mean)
				max.lik.rate <- ovp$MLRate
				max.lik <- ovp$Max.lik
				lik1 <- likRatePhylo(rateData, rep(1,length(rateData$Vmat)), common.mean=common.mean)$ll
				D <- 2 * (max.lik - lik1)
	
	
				mlparams <- likRatePhylo(rateData, rate=max.lik.rate, common.mean=common.mean)
	
	
				if (common.mean==TRUE) { mu <- mlparams$mu } else {
					
					mu <- rep(NA, length(mlparams$mu))
					mu[1] <- mlparams$mu[1]
					for (i in 2:length(mlparams$mu)) { mu[i] <- mlparams$mu[1] + mlparams$mu[i] }
				}
		
														
		
				
				if (length(op) == length(which(op==FALSE)))  { k2 <- length(op) - 1 
						} else { k2 <- length(which(op==FALSE)) }
				
				if(length(op)!=length(which(op==FALSE))) {
					if(common.mean==TRUE) {k <- 2 + length(which(op==FALSE))
							} else { k <- (length(which(op==FALSE)) +1 + length(op)) }
							
							} else {
					
					if(common.mean==TRUE) {k <- 1 + length(which(op==FALSE))
							} else { k <- (length(which(op==FALSE)) + length(op)) }
							}
								
				pval <- 1- pchisq(D, k2)
				n <- length(rateData$y)
				
				CIs.rate <- RatePhylo.allCI(rateData, max.lik.rate, fixed=rep("FALSE", length(rateData$Vmat)), common.mean=common.mean)
				
				aic <- -2 * max.lik + 2 * k
				aicc <- -2 * max.lik + 2 * k + ((2*k*(k+1))/(n-k-1))
				
				if(common.mean==TRUE) {
						aic.rate1 <- -2 * lik1 + 2 * 2
						aicc.rate1 <- -2 * lik1 + 2 * 2 + ((2*2*(2+1))/(n-2-1))
					} else {
						aic.rate1 <- -2 * lik1 + 2 * (1 + length(op))
						aicc.rate1 <- -2 * lik1 + 2 * (1 + length(op)) + ((2*(1 + length(op))*((1 + length(op))+1))/(n-(1 + length(op))-1))
						}
	
				if(pretty == TRUE) {
					cat("____________________________\n")
					cat("Maximum likelihood estimation: rates:\n\n")
					cat("Brownian variance (rate): ", mlparams$s2, "\n")
					cat("ML relative rate estimates:          ", max.lik.rate, "\n")
					cat("Lower confidence intervals for rates:", CIs.rate[,1], "\n")
					cat("Upper confidence intervals for rates:", CIs.rate[,2], "\n")
					cat("Number of parameters: ", k, "\n")
					cat("Maximised log likelihood: ", max.lik, "\n")
					cat("Log likelihood (single rate): ", lik1, "\n")
					cat("LR statistic (ML rates vs single rate):", D)
					cat("  P = ", pval )
					cat(" df = ", k2, " \n")
					cat("  AIC = ", aic, "  \n")
					cat("  AICc = ", aicc, "  \n")
					cat("  Rate 1 AIC = ", aic.rate1, "  \n")
					cat("  Rate 1 AICc = ", aicc.rate1, "  \n")
					cat("____________________________\n")
					cat("Maximum likelihood estimation: means:\n\n")
					cat("ML mean estimates: ", mu, "\n")
					cat("____________________________\n")
					}
				if(pretty == FALSE) {
					max.lik.rate.list <- list(BRvar=mlparams$s2, MLRate = max.lik.rate, LCI = CIs.rate[,1], UCI = CIs.rate[,2], nParam = k, Max.lik = max.lik, Lik1 = lik1, LR = D, P = pval, df = k2, AIC = aic, AICc=aicc, AIC.rate1=aic.rate1, AICc.rate1=aicc.rate1, means=mu)
					return(max.lik.rate.list) }	
	}

