\name{plotROCwithCI}
\alias{plotROCwithCI}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Plot ROC curve with confidence limits
}
\description{
Calculates and plots an annotated ROC curve with confidence limits
}
\usage{
plotROCwithCI(gs, fit3 = fit3, coef, main, pAUCthreshold = 0.2, add = FALSE,
 \dots)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{gs}{a \code{\link[GSEABase]{GeneSet}} object or character vector containing genes-of-intereset
 annotated by Entrez Gene IDs.}
  \item{fit3}{a \code{\link[limma]{eBayes}} object contain comparisons of interest and associated ratios}
  \item{coef}{a \code{numeric} or \code{character} specifying which contrast in \code{fit3} object should be used to
 calculate enrichment.}
  \item{main}{plot title}
  \item{pAUCthreshold}{define cutoff of partial area-under-curve calculation.}
  \item{add}{if set to \code{TRUE}, the ROC curve is palced onto the current plot.}
  \item{\dots}{passes additional functions to \code{plot()} function.}
}
\details{
This function plots an receiver operating characteristic (ROC) curve.

This is a visualization of relative enrichment of gene expression of a given list of genes, specified by \code{gs}, in
 an overall gene expression profile, specified by \code{fit3} and \code{coef}.

The area-under-curve (AUC) and partial area-under-curve (pAUC) are calculated to provide a quantitative value of overall
 enrichment.

This function is a wrapper of functionality provided by and described in \code{\link[ROC]{rocdemo.sca}}.
}
\value{
This function returns a list with two slots:\cr
\enumerate{
 \item \code{AUC} - calculated area under curve
 \item \code{pAUC1} - calculated area under curve
}
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
 \code{\link{plotROC}}
 \code{\link{ROC.compare2gs}}
 \code{\link{ROC.compare2coef}}
 \code{\link[ROC]{rocdemo.sca}}
}
\examples{
data(eset) # load ExpressionSet object 
chipAnnotation <- annotation(eset) # set chip annotation variable

# find enrichment of GO category
library(org.Hs.eg.db)
gs.EG <- unlist(mget("GO:0042063", org.Hs.egGO2ALLEGS)) # retrieve all gliogenesis gene IDs

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

\dontrun{plotROCwithCI(gs.EG, fit3, coef = 1, main = "Enrichment of gliogenesis genes")}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
