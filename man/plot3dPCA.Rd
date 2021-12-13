\name{plot3dPCA}
\alias{plot3dPCA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot 2d and 3d principle component analyses of expression set data
}
\description{
Calculate principle components of each sample in expression data set and makes PCA plots
 using 2d and 3d scatter plots.
}
\usage{
plot3dPCA(eset, groups = NULL, groupnames = NULL, addtext = NULL,
 x.coord = NULL, y.coord = NULL, screeplot = FALSE, squarepca = FALSE,
 pch = NULL, col = NULL, pcs = c(1, 2, 3), legend = TRUE, \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
\item{eset}{either an \code{ExpressionSet} or \code{matrix} of microarray or equivalent data}
\item{groups}{a character \code{factor} containing group assignments for each sample, levels determines
 coloring of groups}
\item{groupnames}{
    a character \code{vector} contain custom group names. Default is \code{NULL}, in which case the
    \code{eset$phenotype} will be used.}
\item{addtext}{
    A character \code{vector} of additional text to be placed close to
    the plotting symbol for each sample. This is helpful if there
    are a lot of samples for identifying e.g., outliers.}
\item{x.coord}{Pass an x-coordinate if automatic legend placement fails}
\item{y.coord}{Pass a y-coordinate if automatic legend placement fails.}
\item{screeplot}{
	Boolean: Plot a \code{\link[stats]{screeplot}} instead of a PCA plot? Defaults to \code{FALSE}.}
\item{squarepca}{Should the y-axis of the PCA plot be made comparable
    to the x-axis? This may aid in interpretation of the PCA
    plot. Defaults to \code{FALSE}.}
\item{pch}{A numeric \code{vector} indicating what plotting symbols to
    use. Default is \code{NULL}, in which case default plotting symbols
    will be used. Note that this argument will override the 'groups'
    argument.}
\item{col}{A numeric or character \code{vector} indicating what color(s) to
    use for the plotting symbols. Default is \code{NULL} in which case
    default colors will be used. Note that this argument will override
    the 'groups' argument.}
\item{pcs}{A character \code{vector} of length two (or three if plot3d
    is \code{TRUE}), indicating which principal components to plot. Defaults to the first two
 principal components.}
\item{legend}{Boolean. Should a legend be added to the plot? Defaults
    to \code{TRUE}.}
\item{\dots}{Further arguments to be passed to \code{plot}. See the help page
    for \code{plot} for further information.}
}
\details{
This function plots both a 3d and 2d PCA of any \code{ExpressionSet} or \code{matrix} object.

If an \code{ExpressionSet} is provided, the function will query the \code{phenotype} data slot of the pData(eset), i.e.
 pData(eset)$phenotype. This will be used to identify grouping within the dataset and coloring of individual samples.

Unlike \code{\link[affycoretools]{plotPCA}}, this function provides support for 3-dimensional PCA plots using the
 \code{\link[rgl]{plot3d}} function. 
}
\value{
 This function returns nothing. It is called only for the side effect of
 producing a PCA plot or screeplot.
}
\references{
 The code for the 2d PCA plot was adapted from \code{\link[affycoretools]{plotPCA}}.
}
\author{
Fraser Sim <fjsim@buffalo.edu>
}
%\note{
%%  ~~further notes~~
%}
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
 \code{\link{customTree}} \cr
 \code{\link[affycoretools]{plotPCA}}
}
\examples{
data(eset)
plot3dPCA(eset)
plot3dPCA(eset, addtext = sampleNames(eset))
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
