\name{pic.motmot}
\alias{pic.motmot}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Phylogenetically independent contrasts (internal)}
\description{
Calculates phylogenetically independent contrasts.  }
\usage{
pic.motmot(x, phy)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
A matrix of trait values with taxon names as rownames.
}
  \item{phy}{
An object of class "phylo" (see ape package).
}

}
\details{
This is a minor tweak to the pic function that extracts the expected variance at the root of the tree (required for likelihood calculations). Here it is only used as an internal function. }
\value{
  \item{contr }{A matrix with two columns containing raw contrasts in the first column and their expected variances in the second column.}
  \item{root.v }{Expected variances of branches either side of root}
  \item{V }{Expected variance at the root }
  \item{ans }{Additional output from pic (mostly not used).}

}
\references{
Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.
}
\author{
Gavin Thomas, Rob Freckleton, Emmanuel Paradis}


%% ~Make other sections like Warning with \section{Warning }{....} ~


\examples{

# Data and phylogeny
data(anolis.tree)
data(anolis.data)

# anolis.data is not matrix and contains missing data so put together matrix of # relevant traits (here female and male snout vent lengths) and remove species 
# with missing data from the matrix and phylogeny
anolisSVL <- data.matrix(anolis.data)[,5]
anolisSVL <- as.matrix(log(anolisSVL))

tree <- drop.tip(anolis.tree, names(attr(na.omit(anolisSVL), "na.action")))
anolisSVL <- na.omit(anolisSVL)

# contrasts
anolis.contr <- pic.motmot(x=anolisSVL, phy=tree)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
