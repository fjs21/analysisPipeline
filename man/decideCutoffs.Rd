\name{decideCutoffs}
\alias{decideCutoffs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Visualize the extent of differential gene expression
}
\description{
This function plots a series of horizontal bar charts to assist in the selection of differential gene expression criteria.
}
\usage{
decideCutoffs(fit3, fc = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{fit3}{an \code{\link[limma]{MArrayLM-class}} object produced by eBayes (limma).}
  \item{fc}{a numeric value specifying the minimum fold change cutoff that will used to identify differentially expressed genes. (not log fold change)}
}
\details{
This function plots a series of bar charts to enable the user to visualize the number of significantly regulated probes using a variety of statistical cutoffs.\cr
In addition, the number of probes above or below a specific fold change is shown.
}
\value{
This function returns nothing. It is called only for the side effect of
 producing decideCutoffs plot.
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
%
%% ~Make other sections like Warning with \section{Warning }{....} ~
\seealso{
\code{\link{parseContrasts}}, \code{\link[limma]{eBayes}}
}
\examples{
# load data
data(eset)

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

decideCutoffs(fit3)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
