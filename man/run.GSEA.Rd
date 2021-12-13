\name{run.GSEA}
\alias{run.GSEA}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
run Gene Set Enrichment 
}
\description{
A single line command to perform GSEA on given expression set and gene set collection objects.
}
\usage{
run.GSEA(nsF, gsc, refSample,
 file = "GSEA.pdf", blnPrompt = FALSE,
 range = c(25, 500), coef = NULL, skip = "NO", pval = 0.05, adjmethod = "BH")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{nsF}{an ExpressionSet object. (Must be filtered to remove duplicately annotated probes.
 featureNames must match the identifiers in the gene set collection object.)}
  \item{gsc}{a Gene Set collection object.}
  \item{refSample}{a \code{character} string specifying which samples to use as control for ploting gene set enrichment.}
  \item{file}{a character file name to store the PDF output.}
  \item{blnPrompt}{a \code{boolean} specifying whether the code waits for the user to proceed following plot commands.}
  \item{range}{specifies the minimum and maximum size of valid gene set to be used in GSEA.}
  \item{coef}{specifies the coeficient of interest in the contrast.matrix object.}
  \item{skip}{not sure.}
  \item{pval}{p-value cutoff for PGSEA/limma analysis (passed to topTable).}
  \item{adjmethod}{adjust method for PGSEA/limma analysis (passed to topTable).}
}
\details{
This function performs parametric gene set enrichment analysis (GSEA) using the \link[PGSEA]{PGSEA} package.

This function looks for enrichment of individual GeneSets within a given \code{\link[GSEABase]{GeneSetCollection}}
 object. Enrichment is calculated against a \code{ExpressionSet} object which is filtered to remove
 duplicated probes to the same Entrez Gene ID (see \code{\link[genefilter]{nsFilter}}). Using the global variable
 design which contains a \code{\link[limma]{limma}} design matrix. Limma is then performed to identify enrichent
 GeneSets. Significantly regulated GeneSets are returned and a heatmap of relative enrichment of significant genesets is
 generated and saved as a PDF file.

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
A \code{data.frame} of significantly enriched or depleted pathways.
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
 \code{\link{runAutomatedGSEA}}
}
\examples{
# load expression data
data(eset)
chipAnnotation <- annotation(eset)

# setup nsF object
library(genefilter)
nsF <- nsFilter(eset, var.filter = FALSE)$eset
probes <- featureNames(nsF)
EG <- getGeneIDFromProbeset(probes)

# use limma to define expression ratios
# refer to limma documentation for details on linear model design.
library(limma)
design = model.matrix(~-1 + eset$phenotype)
colnames(design) <- levels(eset$phenotype)
fit <- lmFit(eset, design)
# contrast.matrix may be set up using multiple contrasts
contrast.matrix <- makeContrasts(A2B5 - uns, levels = design)
fit2 <- contrasts.fit(fit, contrast.matrix)
fit3 <- eBayes(fit2)

# load GeneSet
data(KEGGGOgsc)

# KEGG uses Entrez Gene IDs as identifiers in each gene set
featureNames(nsF) <- EG

# run PGSEA 
result <- run.GSEA(nsF, KEGG.mcs, refSample = "uns",
 file = "GSEA.pdf", blnPrompt = FALSE, range = c(25, 500),
 coef = 1, skip = "YES", pval = 0.05, adjmethod = "none")  
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
