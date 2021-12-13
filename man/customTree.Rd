\name{customTree}
\alias{customTree}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate and plot colored custom hierarchical clustered dendrogram
}
\description{
Calculate and plot colored hierarchical clustering of samples from ExpressionSet object or matrix.
}
\usage{
customTree(eset, hc_labels, groups, title = "Hierarchical Clustering")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{
either an \code{ExpressionSet} or \code{matrix} of microarray or equivalent data.
}
  \item{hc_labels}{
a character vector of labels for each sample in eset.\cr
 defaults to \code{sampleNames(eset)} for \code{ExpressionSet} data\cr
 defaults to \code{colnames(eset)} for \code{matrix} data
}
  \item{groups}{
a character factor containing group assignments for each sample, levels determines coloring of groups.\cr
 defaults to pData(eset)$phenotype for \code{ExpressionSet} data
}
  \item{title}{
a descriptive plot title
}
}
\details{
Plots a hierarchical clustered dendrogram using euclidean distance between samples in
 \code{ExpressionSet} or \code{matrix} of expression data.

Leafs and common branches are colored according to groupings to allow easily visualization
 of sample clustering.
}
\value{
This function returns nothing. It is called only for the side effect of
 producing the dendrogram.
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
  \code{\link{plot3dPCA}}
}
\examples{
data(eset)
customTree(eset)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
% \keyword{ ~kwd1 }
% \keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
