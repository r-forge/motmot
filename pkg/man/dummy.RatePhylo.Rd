\name{dummy.RatePhylo}
\alias{dummy.RatePhylo}
\title{
Trait simulation
}
\description{
Simulates trait data according to a transformed phylogeny. Also allows different means.}
\usage{
dummy.RatePhylo(rateData, rate = NULL, group.means = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{rateData}{
an object of class "rateData"}
  \item{rate}{
a vector of relative rate parameters. The length of the vector is equal to the number of rates 
being estimated. If \code{rate=NULL} then rates are equal.}
  \item{group.means}{
a vector of the relative difference in means between rate categories, expressed as a scalar applied 
to the expected standard deviation (see Ricklefs 2006)}
}
\details{
}
\value{

\item{ydum }{A vector of trait values.}

}
\references{
Ricklefs RE. 2006. Time, species, and the generation of trait variation in clades. Systematic 
Biology 55, 151-159.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments
and island effects. Evolution 63, 2017-2030.}
\author{
Gavin Thomas, Rob Freckleton}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
## Read in phylogeny and data from Thomas et al. (2009)
data(anolis.tree)
data(anolis.data)

## Convert data to class rateData with a rateMatrix object as input
anolis.rateMatrix <- as.rateMatrix(phy=anolis.tree, x=anolis.data$geo_ecomorph, data=anolis.data)

anolis.rateData <- as.rateData(y=anolis.data$Female_SVL, x=anolis.data$geo_ecomorph, 
rateMatrix = anolis.rateMatrix, phy=NULL, data=anolis.data) 

# This is a function for generating dummy data. Traits can be simulated for a given phylogeny and 
# group definition. Both rates and means of the groups can be simulated. Means for different groups
# are generated in relation to the number of expected standard deviations.

dummy.RatePhylo(rateData=anolis.rateData, rate=c(1,1,1,1), group.means=c(0,0,0,0))

dummy.RatePhylo(rateData=anolis.rateData, rate=c(1,2,2,1), group.means=c(0,0,0,0))

dummy.RatePhylo(rateData=anolis.rateData, rate=c(1,1,1,1), group.means=c(0,1,2,3))

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
