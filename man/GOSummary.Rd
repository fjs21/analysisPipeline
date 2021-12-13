\name{GOSummary}
\alias{GOSummary}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performs GO summary on a gene list
}
\description{
This function categorizes differentially expressed genes in common subtypes according to GO annotations.
}
\usage{
GOSummary(geneTable)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{geneTable}{a \code{data.frame} with probe IDs in the named \code{ID} column.}
}
\details{
Differentially expressed genes are first categorized into ligands, receptors, transcription factors, extra-
cellular matrix (ECM), and enzyme genes. The \code{geneTable} is then re-constructed and returned.
}
\value{
a \code{data.frame} with addition of first columned name 'Type' which refers to gene type as either 'ligand',
'receptor', 'transcription factor', 'ECM', 'enzyme' or 'other.
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
%
\seealso{
 \code{\link{runAutomatedGOSummary}}
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# load contrast.details
data(contrast.details)

GOSummary(contrast.details[['A2B5 - uns']]$UPau)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
