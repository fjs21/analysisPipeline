\name{calculateGOHyperG}
\alias{calculateGOHyperG}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perform GO hypergeometric testing
}
\description{
A wrapper function to perform topGO analysis on a given contrast in \code{contrast.details}.
}
\usage{
calculateGOHyperG(contrast, eset, contrast.details)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{contrast}{a \code{character} variable specifying which contrast to analyze.}
  \item{eset}{an \code{ExpressionSet} object of expression data to be analyzed.}
  \item{contrast.details}{a \code{list} object generated by \code{\link{parseContrasts}}.}
}
\details{
Performs \code{\link{topGO}} based analysis of up- and down-regulated genes in the 
specified contrast. \cr The result is written to file, named using the contrast and ending
 'UP_topGO.csv' and 'DOWN_topGO.csv'. 
}
\value{
A list containing results of hypergeometric tests using the topGO function.
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
\code{\link{topGO}}, \code{\link{runAutomatedHyperG}}
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# load contrast.details
data(contrast.details)

\dontrun{contrast.HyperG <- calculateGOHyperG('A2B5 - uns', eset, contrast.details)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
