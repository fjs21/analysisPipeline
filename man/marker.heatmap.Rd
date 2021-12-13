\name{marker.heatmap}
\alias{marker.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plots a marker gene heatmap
}
\description{
This function plots a heatmap of a specific set of marker genes which have been classified according to cell type or pathway involvement.
}
\usage{
marker.heatmap(eset, marker.probes = get("marker.probes", envir = environment(marker.heatmap)),
 clustered = FALSE, title = "Heatmap - Marker Genes",
 groups = NULL, maxLabels = 50, sampleLabels = NULL, geneLabels = TRUE,
 relativeTo = "all", fnHeatmapCol = heatmapCol, lim = NULL)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eset}{an \code{ExpressionSet} object or \code{matrix} containing gene expression data to be plotted.}
  \item{marker.probes}{a \code{character} vector of probe IDs to be plotted on heatmap}
  \item{clustered}{a \code{boolean} specifying whether to cluster samples and genes according to gene expression profiles.}
  \item{title}{the title of the heatmap}
  \item{groups}{a \code{factor} determining the grouping of samples. Used to order samples and to generate a color bar above the heatmap for sample annotation.}
  \item{maxLabels}{a numeric value specifying maximum number of genes, above which no labels will be plotted}
  \item{sampleLabels}{a character vector specifying labels for samples, if not specified sampleNames will be taken from eset$phenotype}
  \item{geneLabels}{a \code{boolean}, if TRUE annotated NCBI Gene Symbols will be used, otherwise rownames of data will be used.}
  \item{relativeTo}{a \code{numeric} vector specifying which samples should be used to normalize gene expression. If not specified, each gene is plotted following normalization against the median expression of all samples. If set to 'none', then no normalization of data is performed.}
  \item{fnHeatmapCol}{a \code{function} to compute the mapping of expression values to colors.\cr
 The default is \code{\link{heatmapCol}} which calculates the upper and lower limits of the data following normalization, and generates a heatmap color bar which is centered around 0, and equally emphasizes gene expression above and below 0.\cr
 The built-in function \code{\link{heatmapCol_stretch}} centers the color bar around 0, and sets the extreme color to the maximal regulation. This is usefull when normalizing data using \code{relativeTo}.}
  \item{lim}{is a numeric value passed to \code{fnHeatmap} which can be used to set manual limits to the color bar.}
}
\details{
This function will plot a heatmap of a given list of probes as so-called \code{marker.probes}. These are typically
 probes whose gene expression is known or expected in a given experiment and can be used to calibrate or assess the
 quality of the underlying gene expression data and experiment.

The global variable \code{markerGenes} should be set as a data.frame prior to running this command. This data.frame
 should contain NCBI Gene Symbols as a \code{character} vector in the \bold{SYMBOL} column, and a \code{numeric}
 vector indicating grouping an position in the heatmap in the \bold{TblOrder} column. This is typically generated in
 an external program like MS Excel and imported into R using \code{read.csv} or similar.

To cluster data by sample and gene expression set the \code{clustered} option to \code{TRUE}.

See \code{\link{custom.heatmap}} for more details on other options.
}
\value{
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
\code{\link{custom.heatmap}}
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# define markerGenes dataframe, 2 columns, SYMBOL = Official Gene Symbosl, TblOrder = Order of groups
markerGenes <- data.frame(cbind(
 c("PDGFRA","CSPG4","ST8SIA1","GFAP","AQP4","S100B"),
 c(1,1,1,2,2,2)))
colnames(markerGenes) <- c("SYMBOL", "TblOrder")
markerGenes

# find annotated probes in data set
library(hgu95av2.db)
marker.probes <- getProbesetFromSymbol(markerGenes$SYMBOL)

# plot unclustered
marker.heatmap(eset)

# plot clustered
marker.heatmap(eset, clustered = TRUE)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
