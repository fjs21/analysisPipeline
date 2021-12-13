\name{bioCarta.heatmap}
\alias{bioCarta.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
plots a heatmap for genes represent in a given bioCarta pathway
}
\description{
This function plots a heatmap for genes represent in a given bioCarta pathway
}
\usage{
bioCarta.heatmap(pwNumber, nsF, species = "Hs")
}
\arguments{
  \item{pwNumber}{pathway number}
  \item{nsF}{expression set object used to retrieve data}
  \item{species}{default = "Hs", human data. Mouse ratios will also work using "Mm"}
}
\details{
This function finds the corresponding probes for a given bioCarta pathway present 
within an ExpressionSet object and generates a heatmap of expression data for 
those probes.

This is a useful way of presenting data from a given pathway identified by, for 
example gene set enrichment analysis (GSEA).

The global variable \code{chipAnnotation}, if not specified, is taken from
 \code{annotation(nsF)}.
The global variable \code{species} allows plotting of human bioCarta-based 
pathways from data generated in other species.

Currently only human "Hs" and mouse "Mm" are supported.

Following annotation of probes to a given pathway, this function calls
 \code{\link{custom.heatmap}}.

To do:
Add functionality for rat and bovine species datasets.
}
\value{
From \code{\link{custom.heatmap}}:
A numeric \code{matrix} of gene expression of selected probes as plotted in heatmap. 
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
 \code{\link{bioCarta.heatmap}}, \code{\link{custom.heatmap}},
 \code{\link{runAutomatedGSEA}}
}
\examples{
 # load example ExpressionSet
 data(eset)
 # filter to remove duplicate probes annotated to same gene on
 # basis of IQR (default)
 library(genefilter)
 nsF <- nsFilter(eset, var.filter = FALSE)$eset

 bioCarta.heatmap(10, nsF) # plots heatmap of mtor signaling pathway

 data(bioCarta)
 print(PathwayNames) # lists all available bioCarta pathways
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
