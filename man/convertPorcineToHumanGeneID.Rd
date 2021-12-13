\name{convertPorcineToHumanGeneID}
\alias{convertPorcineToHumanGeneID}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
convert pig Entrez GeneID to Human Entrez GeneID
}
\description{
This function converts a vector of pig GeneIDs to their human homologs.
}
\usage{
convertPorcineToHumanGeneID(SsEG, verbose = F)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{SsEG}{a vector of pig Entrez Gene IDs}
  \item{verbose}{set to true for debugging purposes}
}
\details{
This function uses three cross-species databases to find human
 homologs of pig genes. \cr

1. Homologene base - uses \code{\link[annotationTools]{getHOMOLOG}} \cr
2. Ensembl biomaRt based - uses \code{\link[biomaRt]{getBM}} \cr

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
# a vector of pig GeneIDs
SsEG = c('407607','397029')
convertPorcineToHumanGeneID(SsEG)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
