\name{fairProportions}
\alias{fairProportions}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Calculate fair proportions phylogenetic diversity metric}
\description{
Calculate fair proportions phylogenetic diversity metric}
\usage{
fairProportions(phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{phy}{
An object of class "phylo" (see ape package).
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{

 \item{fpTips }{Fair proportions metric for all tips in phylogeny}

}
\references{
Redding, D.W. and Mooers, A.O. (2006). Incorporating evolutionary measures into conservation prioritisation. Conservation Biology, 20, 1670-1678.

Isaac, N.J.B., Turvey, S.T., Collen, B., Waterman, C. and Baillie, J.E.M. (2007). Mammals on the EDGE: conservation priorities based on threat and phylogeny. PLoS ONE, 2, e296.}
\author{
Gavin Thomas
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

data(anolis.tree)

fp <- fairProportions(anolis.tree)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
