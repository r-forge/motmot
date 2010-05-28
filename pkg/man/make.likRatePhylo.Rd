\name{make.likRatePhylo}
\alias{make.likRatePhylo}
\title{
Internal function}
\description{
Internal function. Constructor to allow fixing of rate parameters. }
\usage{
make.likRatePhylo(rateData, fixed, common.mean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rateData}{
an object of class "rateData"}
  \item{fixed}{
A vector stating whether each parameter should be allowed to vary (either \code{FALSE} which 
results in a start value of 1, or a numeric start value) or should be fixed (\code{TRUE}).  }
  \item{common.mean}{
a logical specififying whether each rate category should have its own mean (\code{common.mean=FALSE}) 
or all categories should have the same mean (\code{common.mean=FALSE}). See Thomas et al. (2009) 
for a discussion on the impact of assumptions about mean on rate estimates.  }
}
\details{
}
\value{
Returns a function to be passed to \code{optim.likRatePhylo}
}
\references{
}
\author{
Gavin Thomas}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
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
			
			x <- make.anc(y, x)
			
			if(common.mean==FALSE) {
				x <- x} else { x <- rep(1, length(x[,1]))}
						
			logDetV <- logDetFun(retMat)
			iV <- solve(retMat)
			xVix <- crossprod(x, iV \%*\% x)
			xViy <- crossprod(x, iV \%*\% y)
			mu <- solve(xVix) \%*\% xViy 

			e <- y - x \%*\% mu
			s2 <- crossprod(e, iV \%*\% e)
			n <- length(y) 
			phylo.var <- ( s2 / (n - k) )
		
			n <- length(y)
			ll <- -n / 2.0 * log( 2 * pi) - n / 2.0 * log(phylo.var) - logDetV / 2.0 - (n - k)/2.0
			
			ypred <- x\%*\%mu
			lik.RatePhylo <- ( list(ll = ll, mu = mu, phylo.var = phylo.var) )
			return(-1 * lik.RatePhylo$ll)			
	}}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
