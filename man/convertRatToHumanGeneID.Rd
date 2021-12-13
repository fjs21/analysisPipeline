\name{convertRatToHumanGeneID}
\alias{convertRatToHumanGeneID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert rat Entrez GeneID to Human Entrez GeneID
}
\description{
This function converts a vector of rat GeneIDs to their human homologs.
}
\usage{
convertRatToHumanGeneID(RnEG, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{RnEG}{a vector of rat Entrez Gene IDs}
  \item{verbose}{set to true for debugging purposes}
}
\details{
This function uses three cross-species databases to find human
 homologs of rat genes. \cr

1. Homologene base - uses \code{\link[annotationTools]{getHOMOLOG}} \cr
2. Ensembl biomaRt based - uses \code{\link[biomaRt]{getBM}} \cr
3. Bioconductor homology package - \code{\link[hom.Rn.inp.db]{hom.Rn.inp.db}}

Unannotated genes are returned with a GeneID value = 0.

The order and length of the original vector is maintained to permit easy replacment
 of gene indefiers across species.
}
\value{
A vector of human Entrez Gene IDs.
}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Fraser Sim <fjsim@buffalo.edu>
}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{convertHumanToBovineGeneID}},
 \code{\link{convertBovineToHumanGeneID}},
 \code{\link{convertBovineToHumanGeneID}},
 \code{\link{convertHumanToMouseGeneID}},
 \code{\link{convertHumanToRatGeneID}},
}
\examples{
# a vector of rat GeneIDs
RnEG = c('25267','24547')
convertMouseToHumanGeneID(RnEG)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
