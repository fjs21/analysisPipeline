\name{runAutomatedGSEA}
\alias{runAutomatedGSEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
perfomed automated Gene Set Enrichment Analysis
}
\description{
This function performs automated Gene Set Enrichment Analysis (GSEA) on selected contrasts against bioCarta,
Broad C2, Broad C3, GO and KEGG databases.
}
\usage{
runAutomatedGSEA(selCoef, refSample, GSEAfolder = NULL, blnPrompt = FALSE, species = "Hs",
 pcutoff = 0.05, adjmethod = "fdr")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{selCoef}{a \code{character} string specifying the contrast of interest.}
  \item{refSample}{a \code{character} string specifying which samples to use as control for ploting gene set enrichment.}
  \item{GSEAfolder}{a \code{character} string specifying the folder name in which to store the results of GSEA.}
  \item{blnPrompt}{a \code{boolean} specifying whether the code waits for the user to proceed following plot commands.}
  \item{species}{a \code{character} variable specifying the species of the expression set data.\cr
   Only 'Hs' and 'Mm' are currently supported.}
  \item{pcutoff}{p-value cutoff.}
  \item{adjmethod}{a \code{character} variable specifying the adjmethod used by limma for calculation of enriched gene sets.\cr
   This is passed to \code{\link[limma]{topTable}}. The default is false discovery rate, 'fdr'.}
}
\details{
This function automates the process of performing parametric gene set enrichment (PGSEA) on a number of
 bioinformatics databases.

Using \code{\link[PGSEA]{PGSEA}} as the basis for gene set enrichment and \code{\link[limma]{limma}} to calculate
 significance based on the experimental design specified in \code{design}.

See \code{\link{run.GSEA}} for details on individual GSEA analysis.

The following GeneSetCollections are analyzed in turn to identify enriched/depleted GeneSets:
\enumerate{
 \item{} bioCarta - data parsed from \url{www.bioCarta.com}
 \item{} broadC2 - to enable BroadC2 - the user must first download the C2.all.v3.0.symbols.gmt file from
 Broad mSigDb and enter the file location in the global variable \code{GSEA.C2}. Once this variable is set,
 subsequent calls to \code{runAutomatedGSEA} will utilize the BroadC2 database. See vignette for details.
 \item{} broadC3 - to enable BroadC3 - the user mustfirst download the C3.all.v3.0.symbols.gmt file from
 Broad mSigDb and enter the file location in the global variable \code{GSEA.C3}. 
 \item{} GO
 \item{} KEGG
}

GSEA.C2 - \url{http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/3.0/c2.all.v3.0.symbols.gmt}\cr
GSEA.C3 - \url{http://www.broadinstitute.org/gsea/msigdb/download_file.jsp?filePath=/resources/msigdb/3.0/c3.all.v3.0.symbols.gmt}

\code{runAutomatedGSEA} currently supports human and mouse data only.

To enable identification of genesets enriched or depleted among multiple contrasts, the user can specify a
 \code{desiredResults} matrix prior to running \code{run.GSEA}. This matrix should contain the same number of
 columns as contrasts are present in the global variable \code{contrast.matrix}. Each row conditions a logic
 value to select regulated genesets. \cr
A value of 1 specifies that the geneset must be enriched, 0 = no change, -1 = depleted, NA = any regulation allowed.

For example:
\tabular{cccc}{
 \bold{Contrast} \tab A - B \tab A - C \tab B - C \cr
 Row 1 \tab 1 \tab 1 \tab NA \cr
 Row 2 \tab -1 \tab -1 \tab NA
}
In the above example, run.GSEA will identify only those genesets which were enriched in both A vs. B and A vs. C
 (row 1), and those genesets depleted in both A vs. B and A vs. C.

There is no limitation of rows or combinations of contrasts.
}
\value{
A series of csv and PDF files are generated when running this function.
No values are returned.
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
 \code{\link{runAutomatedHyperG}}, \code{\link{runAutomatedGOSummary}}
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

# run automated GSEA on first contrast - uncomment next line to run
\dontrun{runAutomatedGSEA(selCoef = "A2B5 - uns", refSample = "uns", blnPrompt = TRUE)}

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
