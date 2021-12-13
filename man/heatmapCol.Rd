\name{heatmapCol}
\alias{heatmapCol}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to calculate heatmap color bar
}
\description{
Calculates color bar for \code{\link{custom.heatmap}} and \code{\link{marker.heatmap}} functions.
}
\usage{
heatmapCol(data, col, lim)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{a numeric \code{matrix} containing gene expression values.}
  \item{col}{a color range on which to map gene expression values.}
  \item{lim}{the limit value which will be used to determine mapping to extreme colors.
 This is typically calculated using \code{min(abs(range(data)))}.}
}
\details{
This function is not usually called directly but as a consequence of running \code{\link{custom.heatmap}}.

This function calculates a calculates the colors to plot a gene-expression heatmap.\cr
\code{heatmapCol} maps colors to the minimal extent of the range of data. This provides clear visualization
 of gene expression changes from the median of the data.
}
\value{
A \code{character} vector of colors for use in plotting heatmaps.
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
 \code{\link{heatmapCol_stretch}}, \code{\link{custom.heatmap}}, \code{\link{marker.heatmap}}
}
\examples{
# load data
data(eset)

data <- exprs(eset)[1:50,] # take data from first fifty genes
data.medians = apply(data,1,median)
data <- sweep(data, 1, data.medians)

hmcol = colorRampPalette(brewer.pal(10, "RdYlGn")[10:1])(256) # a typical green-red heatmap color bar
lim = min(abs(range(data))) # calculate limits

heatmapCol(data = data, col = hmcol, lim = lim)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
