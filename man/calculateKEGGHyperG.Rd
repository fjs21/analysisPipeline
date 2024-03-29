\name{calculateKEGGHyperG}
\alias{calculateKEGGHyperG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform KEGG-based hypergeometric testing
}
\description{
A function to perform hypergeometric testing of a gene list on the KEGG pathway database.
}
\usage{
calculateKEGGHyperG(contrast, GeneIDUniverse, contrast.details)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contrast}{a \code{character} variable specifying which contrast to analyze.}
  \item{GeneIDUniverse}{a \code{character} variable containing all NCBI Entrez Gene IDs for the data set.}
  \item{contrast.details}{a \code{list} object generated by \code{\link{parseContrasts}}.}
}
\details{
Performs KEGG-based hypergeometric testing on the up- and down-regulated genes in the specified contrast.

This function looks for over-representation of KEGG pathway constituents using a
 default p-value cut-off of 0.05.

The global variable \code{chipAnnotation} determines annotation of probes to unique NCBI Entrez Genes.
Over-representation is determined against the character vector \code{GeneIDUniverse}.

The result is written to file, named using the contrast and ending
 'UP_topGO.csv' and 'DOWN_topGO.csv'. 
}
\value{
A list containing results of hypergeometric tests using the KEGG pathways.
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
 \code{\link{runAutomatedHyperG}}
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# load contrast.details
data(contrast.details)

# define GeneIDUniversse (all genes on chip)
library(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
GeneIDUniverse = unique(unlist(mget(featureNames(eset),
 eval(parse(text = paste (chipAnnotation,"ENTREZID",sep=""))), ifnotfound = NA )))
GeneIDUniverse = GeneIDUniverse[!is.na(GeneIDUniverse)]

# perform KEGG hypergeometric tests
contrast.HyperG <- calculateKEGGHyperG('A2B5 - uns', GeneIDUniverse, contrast.details)

# return table of significant pathways in up-regulated genes
summary(contrast.HyperG$UP_KEGG)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
