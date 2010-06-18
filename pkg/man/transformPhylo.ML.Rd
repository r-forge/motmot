\name{transformPhylo.ML}
\alias{transformPhylo.ML}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Maximum likelhood for models of trait evoluion.}
\description{
Fits likelihood models for various models of continuous character evolution. Model fitting is based on maximum-likelihood evaluation using phylogenetically independent contrasts. This is exactly equivalent to, but substantially faster than, GLS approaches.
}
\usage{
transformPhylo.ML(y, phy, model = NULL, nodeIDs = NULL, minCladeSize=10, nSplits=10, rateMin = 1e-05, rateMax = 1e+05, deltaMax = 10)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{y}{
A matrix of trait values.
}
  \item{phy}{
An object of class "phylo" (see ape package).
}
  \item{model}{
The model of trait evolution (see details).
}
  \item{nodeIDs}{
Integer - ancestral nodes of clades.
}
  \item{minCladeSize}{
Integer - minimum clade size for inferred rate shift (where model="medusa").
}
  \item{nSplits}{
Integer - number of rate shifts to apply for model="medusa"}
  \item{rateMin}{
Minimum value for relative rate estimate}
  \item{rateMax}{
Maximum value for relative rate estimate}
  \item{deltaMax}{
Maximum value for delta}
}
\details{
This function finds the maximum likelihood parameter values for continuous character evolution. The function returns the maximum-likihood parameter estimates for the following models.

model="bm"- Brownian motion (constant rates random walk)

model="kappa" - fits Pagel's kappa by raising all branch lengths to the power kappa. As kappa approaches zero, trait change becomes focused at branching events. For complete phylogenies, if kappa approaches zero this infers speciational trait change. 

model="lambda" - fits Pagel's lambda to estimate phylogenetic signal by multiplying all internal branches of the tree by lambda, leaving tip branches as their original length (root to tip distances are unchanged);

model="delta" - fits Pagel's delta by raising all node depths to the power delta. If delta <1, trait evolution is concentrated early in the tree whereas if delta >1 trait evolution is concentrated towards the tips. Values of delta above one can be difficult to fit reliably.

model="free" - fits Mooer's et al's free model where each branch has its own rate of trait evolution. This can be a useful exploratory analysis but it is slow due to the number of parameters, particularly for large trees.

model="clade" - fits a model where particular clades are a priori hypothesised to have different rates of trait evolution (see O'Meara et al. 2006; Thomas et al. 2006, 2009). Clades are specified using nodeIDs and are defined as the mrca node.

model="medusa" - fits "clade" models without any a priori assertion of the location of phenotypic diversification rate shifts. It uses the same approach as the runMedusa function in the geiger package (runMedusa tests for shifts in the rate of lineage diversification). The algorithm first fits a constant-rate Brownian model to the data, it then works iteratively through the phylogeny fitting a two-rate model at each node in turn. Each two-rate model is compared to the constant rate model and the best two-rate model is retained. Keeping the location of this rate shift intact, it then repeats the procedure for a three-rate model and so on. The maximum number of rate shifts can be specified a priori using nSplits. Limits can also be applied to the size (species richness) of clades on which to infer new rate shifts using minCladeSize. This can be useful to enable large trees to be handled but should be used cautiously since specifiying a large minimum clade size may result in biologically interesting nested rate shifts being missed. Equally, very small clade sizes may provide poor estimates of rate that may not be informative.

}
\value{
  \item{maximum }{The maximum likelihood parameter estimates}
  \item{objective }{The maximum log-likelihood}
  
  \item{traitMedusaObject }{A list in which the first element contains a matrix summarising the parameter estimates and node ids, log-likelihoods, number of parameters (k), AIC and AICc for the best one-rate model, two-rate model, three rate model and so on. The second element is a sub-list where the first element contains all two-rate models, the second element contains all three-rate models and so on. This can be summarised using traitMedusaSummary.}
}
\references{
Alfaro ME, Santini F, Brock CD, Alamillo H, Dornburg A, Carnevale G, Rabosky D & Harmon LJ. 2009. Nine exceptional radiations plus high turnover explain species diversity in jawed vertebrates. PNAS 106, 13410-13414.

Felsenstein J. 1973. Maximum-likelihood estimation of evolutionary trees from continuous characters. Am. J. Hum. Genet. 25, 471-492.

Felsenstein J. 1985. Phylogenies and the comparative method. American Naturalist 125, 1-15.

Freckleton RP & Jetz W. 2009. Space versus phylogeny: disentangling phylogenetic and spatial signals in comparative data. Proc. Roy. Soc. B 276, 21-30. 

Mooers AO, Vamosi S, & Schluter D. 1999. Using phylogenies to test macroevolutionary models of trait evolution: sexual selection and speciation in Cranes (Gruinae). American Naturalist 154, 249-259.

O'Meara BC, Ane C, Sanderson MJ & Wainwright PC. 2006. Testing for different rates of continuous trait evolution using likelihood. Evolution 60, 922-933

Pagel M. 1997. Inferring evolutionary processes from phylogenies. Zoologica Scripta 26, 331-348.

Pagel M. 1999 Inferring the historical patterns of biological evolution. Nature 401, 877-884.

Thomas GH, Meiri S, & Phillimore AB. 2009. Body size diversification in Anolis: novel environments and island effects. Evolution 63, 2017-2030.
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

# Data and phylogeny
data(anolis.tree)
data(anolis.data)

# Resolve polytomies
tree <- multi2di(anolis.tree)

# anolis.data is not matrix and contains missing data so put together matrix of # relevant traits (here female and male snout vent lengths) and remove species 
# with missing data from the matrix and phylogeny
anolisSVL <- data.matrix(anolis.data)[,c(5,6)]
anolisSVL[,1] <- log(anolisSVL[,1])
anolisSVL[,2] <- log(anolisSVL[,2])

tree <- drop.tip(tree, names(attr(na.omit(anolisSVL), "na.action")))
anolisSVL <- na.omit(anolisSVL)

# Brownian motion model
transformPhylo.ML(anolisSVL, phy=tree, model="bm")

# Pagel's kappa, lambda, and delta
transformPhylo.ML(anolisSVL, phy=tree, model="kappa")
transformPhylo.ML(anolisSVL, phy=tree, model="lambda")
transformPhylo.ML(anolisSVL, phy=tree, model="delta")


# Test for different rates in different clades - here with 2 hypothesised unusual rates compared to the background
# This fits the non-censored model of O'Meara et al. (2006)
transformPhylo.ML(anolisSVL, phy=tree, model="clade", nodeIDs=c(178, 199))

# Identify rate shifts and print and plot results
anolisSVL_MEDUSA <- transformPhylo.ML(anolisSVL, phy=tree, model="medusa", minCladeSize=20, nSplits=3)

anolisSVL_MEDUSA_out <- traitMedusaSummary(anolisSVL_MEDUSA, phy=tree)


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
