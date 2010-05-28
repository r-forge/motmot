\name{make.anc}
\alias{make.anc}
\title{
Internal function
}
\description{
This is an internal function to generate the design matrix required to define different means for 
each hypothesised rate category.}
\usage{
make.anc(y, x, data = NULL, common.mean = FALSE)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
The response variable - typically a continuous trait.}
  \item{x}{
The explanatory (discrete) variable used to define the hypothesised rate categories.}
  \item{data}{
A data frame containing (minimally) the x and y variables as columns with species names as rownames.}
  \item{common.mean}{
a logical specififying whether each rate category should have its own mean 
(\code{common.mean=FALSE}) or all categories should have the same mean (\code{common.mean=FALSE}). 
See Thomas et al. (2009) for a discussion on the impact of assumptions about mean on rate estimates.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

\item{x }{A design matrix}

}
\references{
Thomas GH, Freckleton RP, & Szekely T. 2006. Comparative analyses of the influence of developmental
mode on phenotypic diversification rates in shorebirds. Proceedings of the Royal Society B 273, 
1619-1624.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments
and island effects. Evolution 63, 2017-2030.}
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
function(y, x, data=NULL, common.mean=FALSE) {
	
	if(is.factor(x) == "FALSE"){
		stop("The discrete trait must be a factor")	
	}
	
	m <- model.frame(y ~ x, data)
	x <- model.matrix(y ~ x, m)
	return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
