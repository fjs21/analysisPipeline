\name{parseContrasts}
\alias{parseContrasts}
\title{
Process the results of limma into contrast.details
}
\description{
This function processes the results of a linear model differential gene expression analysis. 
}
\usage{
parseContrasts(fit3, FCcutoff = 2, FDRcutoff = 0.05, adjust = "fdr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit3}{
an \code{\link[limma]{MArrayLM-class}} object produced by eBayes (limma).
}
  \item{FCcutoff}{
a numeric value specifying the minimum fold change cutoff that will used to identify differentially
 expressed genes. (not log fold change)
}
  \item{FDRcutoff}{
a numeric value specifying the p-value cutoff
}
  \item{adjust}{
method used to adjust the p-values for multiple testing. see /link[limma]{topTable} for more details.
}
}
\details{
parseContrasts automates the process of calling topTable on every contrast in a linear model analysis.

This function determines the lists of genes up-regulated and down-regulated in each comparison
 (contrast) and generates summary tables for each as .csv files in the working directory.

The returned \code{contrast.details} can be used for hypergeometric analysis using KEGG and GO with
 \code{\link{runAutomatedHyperG}}.

This function acts like \code{\link[affycoretools]{limma2annaffy}} to generate data tables (in this case as csv
 files) of the significantly regulated genes in each contrast present in \code{fit3}. Unlike \code{\link[affycoretools]{limma2annaffy}}, this
 function supports any platform with a valid annotation package as specified by global variable \code{chipAnnotation}.
}
\value{
Produces a list, with one member per contrast in \code{fit3} object.\cr
Each member contains the following objects:\cr
  \item{topT}{result of \code{\link[limma]{topTable}} using specified settings}
  \item{UP}{a \code{data.frame} containing all up-regulated probes.}
  \item{UPa}{a \code{data.frame} containing all up-regulated probes, probes have been annotated with Symbol, Description, and NCBI Entrez Gene.}
  \item{UPau}{a \code{data.frame} containing all up-regulated probes. Duplicated probes to the same gene have been removed with the most significant probe retained in each case.}
  \item{DOWN}{a \code{data.frame} containing all down-regulated probes.}
  \item{DOWNa}{a \code{data.frame} containing all down-regulated probes, probes have been annotated with Symbol, Description, and NCBI Entrez Gene.}
  \item{DOWNau}{a \code{data.frame} containing all down-regulated probes. Duplicated probes to the same gene have been removed with the most significant probe retained in each case.}
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
 \code{\link{decideCutoffs}}
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# use limma to define differentially expressed transcipts
# refer to limma documentation for details on linear model design.
library(limma)
design = model.matrix(~-1 + eset$phenotype)
colnames(design) <- levels(eset$phenotype)
fit <- lmFit(eset, design)

# contrast.matrix may be set up using multiple contrasts
contrast.matrix <- makeContrasts(A2B5 - uns, levels = design)

fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)

contrast.details <- parseContrasts(fit3) # calculate gene lists using default settings

contrast.details <- parseContrasts(fit3, FCcutoff = 5, FDRcutoff = 0.01, adjust = "fdr") # calculate gene lists using 5 fold change and 1% FDR cutoffs.
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
