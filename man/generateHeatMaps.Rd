\name{generateHeatMaps}
\alias{generateHeatMaps}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to generate heatmaps from each contrast in a contrast.details object
}
\description{
This function plots heatmaps of the top 50 genes for each contrast specified in \code{contrast.details}.
}
\usage{
generateHeatMaps(eset, export = TRUE, contrast.names = NULL, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{export}{a \code{boolean}. If \code{TRUE}, heat maps are saved to PDF file.
 If \code{FALSE}, heatmap are plotted within R.}
  \item{contrast.names}{a \code{character} vector specifying which contrasts to plot.
 If \code{NULL}, all contrasts are plotted.}
  \item{eset}{an \code{ExpressionSet} object containing expresion data to plotted on heatmaps.
}
  \item{\dots}{addition parameters passed onto \code{\link{custom.heatmap}}.
 This can be useful for defining custom groups and color bar schemes.}
}
\details{
This function automates the process of generating heatmaps for each set of up- and down-regulated genes in
 each contrast. The top 50 most significant probes UP and DOWN from each contrast are plotted seperately.

This function relies on a valid \code{contrast.details} object.
}
\value{
This function returns nothing. It is called only for the side effect of
 producing a preformatted heat maps for each contrast in \code{contrast.details}.
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
 \code{\link{custom.heatmap}}, \code{\link{parseContrasts}}
}
\examples{
# load data
data(eset)

# load contrast.details
data(contrast.details)

generateHeatMaps(eset = eset, export = FALSE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
