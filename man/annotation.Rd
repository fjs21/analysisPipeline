\name{annotation}
\alias{getSymbolFromProbeset}
\alias{getGeneIDFromProbeset}
\alias{getDescriptionFromProbeset}
\alias{getProbesetFromSymbol}
\alias{getProbesetFromGeneID}
\title{
simple functions to annotate probe sets and retrieve annotations
}
\description{
A series of functions with single parameters for annotation.
}
\usage{
getSymbolFromProbeset(probes)
getGeneIDFromProbeset(probes)
getDescriptionFromProbeset(probes)
getProbesetFromSymbol(syms)
getProbesetFromGeneID(EG)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{a \code{character} vector of probe sets.}
  \item{syms}{a \code{character} vector of Entrez Gene symbols.}
  \item{EG}{a \code{character} vector of Entrez Gene IDs.}
}
\details{
These functions are designed to simply the annotation of probes to and from Entrez Gene ID and Entrez Gene Symbols.
 They are used throughout the package as simple wrapper functions.

The global variable \code{chipAnnotation} must be set before calling these functions.\cr
This is typically manually set at the begining of the analysis pipeline using \code{chipAnnotation = annotation(eset)}.
}
\value{
\item{getSymbolFromProbeset(probes)}{A named \code{character} vector of gene symbols.
 If no annotation is available 'NA' is returned.}

\item{getGeneIDFromProbeset(probes))}{A named \code{numeric} vector of Entrez Gene IDs.
 If no annotation is available 0 is returned.}

\item{getDescriptionFromProbeset(probes))}{A named \code{character} vector of Entrez gene descriptions.
 If no annotation is available 'NA' is returned.}

\item{getProbesetFromSymbol(syms))}{A named \code{character} vector of probes.
 If no annotation is available 'NA' is returned.}

\item{getProbesetFromGeneID(EG))}{A named \code{character} vector of probes.
 If no annotation is available 'NA' is returned.}
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
%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
\examples{
# load data
data(eset)
chipAnnotation = annotation(eset)

# retrieve annotations
getSymbolFromProbeset(featureNames(eset)[1:20])
getGeneIDFromProbeset(featureNames(eset)[1:20])
getDescriptionFromProbeset(featureNames(eset)[1:20])

# find matching probe sets from SYMBOLS
goi = c("PDGFRA","MOG","MAG")
getProbesetFromSymbol(goi)

# find matching probe sets from GeneIDs
goi = 5156
getProbesetFromGeneID(goi)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
