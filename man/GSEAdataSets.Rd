\name{GSEAdataSets}
\alias{GO.mcs}
\alias{KEGG.mcs}
\alias{InteractionCondition}
\alias{InteractionIDs}
\alias{InteractionMap}
\alias{InteractionType}
\alias{MoleculeIDs}
\alias{MoleculeMap}
\alias{PathwayIDs}
\alias{PathwayMap}
\alias{PathwayNames}
\docType{data}
\title{
Gene Set Collection data data
}
\description{
bioCarta, KEGG and GO data sets for GSEA analysis.
}
%\usage{data(KEGGGOgsc.rda)}
%\format{
%}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
\source{
bioCarta - downloaded xml file from \url{http://www.biocarta.com/} and parsed for analysis

GO.mcs - generated from \code{\link[GO.db]{GO.db}} using \code{\link[PGSEA]{kegg2smc}}

KEGG.mcs - generated from \code{\link[KEGG.db]{KEGG.db}} using \code{\link[PGSEA]{go2smc}}
}
%\references{
%%  ~~ possibly secondary sources and usages ~~
%}
%\examples{
%data(KEGGGOgsc.rda)}
%data(bioCarta.rda)}
%
%}
\keyword{datasets}
