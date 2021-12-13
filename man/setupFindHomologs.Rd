\name{setupFindHomologs}
\alias{setupFindHomologs}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
setup FindHomologs functions
}
\description{
This function loads homologene and connects to biomaRt.
}
\usage{
setupFindHomologs()
}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%}
\details{
This function is required prior to using any of cross-species annotation functions.

In order to utilize the \code{\link[annotationTools]{getHOMOLOG}} function from the annotationTools package
 for cross-species annotations, you will need to download homologene file from NCBI and extract it to a
 specific file location. \cr
(Homologene.data is available here \url{ftp.ncbi.nih.gov/pub/HomoloGene/current/})

Once downloaded, the \code{homologeneFile} global variable should be set to the location of the homologene file
 before executing.
}
\value{
Nothing.
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
 \code{\link{convertHumanToBovineGeneID}},
 \code{\link{convertBovineToHumanGeneID}},
 \code{\link{convertHumanToMouseGeneID}},
 \code{\link{convertMouseToHumanGeneID}},
 \code{\link{convertHumanToRatGeneID}},
 \code{\link{convertRatToHumanGeneID}}
}
\examples{
\donttest{
# run the code below to download the homologene file for cross-species annotation using annotationTools
ftp <- "ftp://ftp.ncbi.nih.gov/pub/HomoloGene/current/homologene.data"
download.file(url = ftp, destfile = "homologene.data")
homologeneFile <- "homologene.data"

# setupFindHomologs
setupFindHomologs()
}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
