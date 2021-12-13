\name{ROC.compare2coef}
\alias{ROC.compare2coef}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ROC curves for two coefficients in same limma analysis and compare enrichment.
}
\description{
Calculates and plots an annotated ROC curve for two defined coefficients. Calculates probability that the two coefficients exhibit different expression profiles for the given gene set.
}
\usage{
ROC.compare2coef(gs, fit3 = fit3, coef1, coef2, main, pAUCthreshold = 0.2, add = FALSE,
 \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gs}{a \code{\link[GSEABase]{GeneSet}} object or character vector containing genes-of-intereset annotated by Entrez Gene IDs.}
  \item{fit3}{a \code{\link[limma]{eBayes}} object contain comparisons of interest and associated ratios}
  \item{coef1}{a \code{numeric} or \code{character} specifying which contrast in \code{fit3} object should be used to calculate enrichment.}
  \item{coef2}{a \code{numeric} or \code{character} specifying which contrast in \code{fit3} object should be used to calculate enrichment.}
  \item{main}{plot title}
  \item{pAUCthreshold}{define cutoff of partial area-under-curve calculation.}
  \item{add}{if set to \code{TRUE}, the ROC curve is palced onto the current plot.}
  \item{\dots}{passes additional functions to \code{plot()} function.}
}
\details{
This function plots receiver operating characteristic (ROC) curve for two gene sets in a given gene expression profile defined via a limma contrast.

This is a visualization of relative enrichment of gene expression of a given list of genes, specified by \code{gs}, in an overall gene expression profile, specified by \code{fit3} and \code{coef}.

The area-under-curve (AUC) and partial area-under-curve (pAUC) are calculated to provide a quantitative value of overall enrichment.

This function is a wrapper of functionality provided by and described in \code{\link[pROC]{plot.roc}}.
}
\value{
This function returns a list with two slots:\cr
\enumerate{
 \item \code{pval} - p-value
 \item \code{AUC1} - calculated area under curve for gene set 1
 \item \code{pAUC1} - calculate partial area under curve for gene set 1
 \item \code{AUC2} - calculated area under curve for gene set 2
 \item \code{pAUC2} - calculate partial area under curve for gene set 2
}
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
 \code{\link{plotROC}}
 \code{\link{plotROCwithCI}}
 \code{\link{ROC.compare2gs}}
 \code{\link[ROC]{rocdemo.sca}}
}
\examples{
# to do
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
