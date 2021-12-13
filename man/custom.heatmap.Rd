\name{custom.heatmap}
\alias{custom.heatmap}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to plot a customizable gene expression heatmap
}
\description{
This function can be used to quickly and easily plot heatmap of specific genes in an
 \code{ExpressionSet} or \code{matrix} object.
}
\usage{
custom.heatmap(heatmap.probes, eset, title = "Heatmap", groups = NULL, maxLabels = 50,
 sampleLabels = NULL, geneLabels = TRUE, relativeTo = "all", fnHeatmapCol = heatmapCol, lim = NULL,
 dendrogram = "row", Colv = FALSE,
 \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{heatmap.probes}{a \code{character} vector of probe IDs to be plotted on heatmap}
  \item{eset}{an \code{ExpressionSet} object or \code{matrix} containing gene expression data to be plotted.}
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
  \item{dendrogram}{passed to \code{heatmap.2()} }
  \item{Colv}{passed to \code{heatmap.2()} }
  \item{\dots}{passes additional functions to \code{heatmap.2()} function.}
}
\details{
This function uses \code{\link[gplots]{heatmap.2}} to plot a gene-expression heatmap. 

The function has been tailored to generate heatmaps appropriately formatted to visualize microarray data
 from a selected list of probes. 

When using an \code{ExpressionSet} object:\cr
The default function uses groupings defined in the \code{pData(eset)$phenotype}. Samples are ordered according to the 
 factor levels in this variable. As such, it is possible to easily order the samples of a complex experiment.\cr
Sample labels are taken from the eset as \code{sampleLabels(eset}.\cr
Gene labels are derived from the eset annotation as \code{annotation(eset)}. 

When using a \code{matrix} object:\cr
If \code{groups = NULL} then no grouping are defined.
The global variable \code{chipAnnotation} must be defined as a valid annotation package, such as "hgu95av2".

By default, expression data is normalized to the median of each probe prior to plotting. This behavior can be 
altered by specifying either \code{relativeTo == "none"} or providing either a logical/boolean variable or
 a vector of selected samples against which to normalize.

By default, the coloring bar is design to maximize visualization of gene expression changes among the majority
 of probes. The \code{lim} variable can be used to reduce the saturation values at which point the colors are either
 ful red or green.

In order to visualize the extreme changes in gene expression, use \code{fnHeatmapCol = \link{heatmapCol_stretch}}.
 This function sets the limits of the color bar to the far extreme of the data set.
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
\code{\link{marker.heatmap}}
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# example 1 - GO pathway
library(hgu95av2.db)
gs.probes <- unlist(mget("GO:0042063", hgu95av2GO2ALLPROBES)) # retrieve all gliogenesis probe sets
custom.heatmap(gs.probes, eset, "Gliogenesis genes")

# example 2 - top 50 genes in contrast
# load contrast.details
data(contrast.details)

gs.probes <- contrast.details[['A2B5 - uns']]$UP$ID[1:50]
custom.heatmap(gs.probes, eset, "Genes up-regulated in A2B5 cells")

# example 3 - use a streched color scale & normalize to 'uns' samples
custom.heatmap(gs.probes, eset, "Genes up-regulated in A2B5 cells", 
 relativeTo = which(eset$phenotype == "uns"),
 fnHeatmapCol = heatmapCol_stretch)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
