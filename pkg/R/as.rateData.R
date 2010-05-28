as.rateData <-
function(y, x, rateMatrix=NULL, phy, data) {
		
		if (is.null(rateMatrix))	{ rateMatrix <- as.rateMatrix(phy=phy, x=x, data=data) } else { rateMatrix <- rateMatrix }
		
		dat <- data.frame(x, y, row.names = rownames(data))
	
		if(sum(sort(rownames(dat)) != sort(rownames(rateMatrix[[1]]))) > 0){warning("Length or names of phenotypic and of phylogenetic data do not match - non-matching taxa will be dropped")}	# check names

		
		# get names of species common to data frame and phylogeny
		sharedSpecies <- intersect(rownames(rateMatrix[[1]]), rownames(dat))
		
		# subset data to only those species in phylogeny and alphabetise
		dat <- dat[match(sharedSpecies ,rownames(dat)),] 
		dat<- dat[sort(rownames(dat), index.return = TRUE)$ix, ]
		
		# index missing data for y
		idx <- which( is.na(dat$y) == FALSE)
	
		# alphabetise and prune rate matrices to remove data with missing species
		Vmat <- vector(mode="list", length = length(rateMatrix))

		for(i in 1:length(rateMatrix)) {
			Vmatrix <- rateMatrix[[i]]
			nms <- rownames(Vmatrix)
			snms <- sort(nms, index.return = TRUE)
			Vmatrix <- Vmatrix[snms$ix, snms$ix]
			Vmatrix <- Vmatrix[idx, idx]
			Vmat[[i]] <- Vmatrix
			}
		
			x <- as.matrix(dat$x)
			y <- dat$y[idx]
			x <- x[idx, ]
			traits <- list(y = y, x = x, Vmat = Vmat)
			class(traits) <- "rateData"
			return(traits)
			}

