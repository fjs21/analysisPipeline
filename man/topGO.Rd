\name{topGO}
\alias{topGO}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
performed topGO analysis and saves result in CSV file
}
\description{
This function performs a topGO analysis and saves the results in a CSV file.
}
\usage{
topGO(probes, eset = eset)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{probes}{a \code{character} vector of probes of interest}
  \item{eset}{an \code{ExpressionSet} object defining the entire gene universe and annotation}
}
\details{
This wrapper function performs topGO analysis using \code{\link[topGO]{runTest}} on a given set of probes in an
 \code{ExpressionSet} object.

The 'elim' algorithm and fisher statistic are used to test for over-representated GO categories among biological
 process terms. The significance cut-off is set at 0.01.

The results of topGO analysis is returned following annotation of terms with significantly regulated genes, as
 both Entrez Gene IDs and Gene Symbols. 

If the expression set is not human, you will need to set mapping to the appropriate species.\cr
E.g. for mouse \code{mapping = 'org.Mm.eg'}, for rat \code{mapping = 'org.Rn.eg'}.
}
\value{
a \code{data.frame} of enriched GO terms.
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
 \code{\link{runAutomatedHyperG}}, \code{\link{calculateGOHyperG}}
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# load contrast.details
data(contrast.details)

\dontrun{topGO(contrast.details[['A2B5 - uns']]$UP$ID, eset)}
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
