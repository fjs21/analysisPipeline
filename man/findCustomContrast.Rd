\name{findCustomContrast}
\alias{findUP}
\alias{findDOWN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
function to process custom lists of UP and DOWN regulated genes in fit3 object
}
\description{
These functions allow the definition of custom contrasts (typically via Venn Diagram analysis)
 for downstream analysis.
}
\usage{
findUP()
findDOWN()
}
\details{
\code{findUP} and \code{findDOWN} are used to define custom contrasts after the initial analysis is complete.
This is useful when perform venn analysis looking for overlap between contrasts.\cr
To use \code{findUP} and \code{findDOWN}, the following global variables must be defined:\cr
\code{contrast.details}, \code{contrast}, \code{contrastRatio}, \code{contrastProbes}

\enumerate{
 \item{} \code{contrast} defines the name of the new contrast.
 \item{} \code{contrastProbes} defines the probes to be considered in the new contrast.
 \item{} \code{contrastRatio} defines the contrast from which data will be taken to generate new table.
 \item{} \code{contrast.details} - the existing contrast.details object generated typically
 following \code{\link{parseContrasts}}
}
}
\value{
\code{contrast.detail} list with addition of new contrast.
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
 \code{\link{parseContrasts}}, ~~~
}
\examples{
# load data
data(eset)
# set chipAnnotation
chipAnnotation = annotation(eset)

# load contrast.details
data(contrast.details)

# find overlapping lists
# select name for new contrast
contrast <- "newContrastName"
# select contrast which will be used to display fold change, significance, etc in tables
contrastRatio <- 1 
# define list of up-regulated probes in new contrast
contrastProbes <- sample(contrast.details[['A2B5 - uns']]$UP$ID,50) # random selection of 50 genes
contrast.details <- findUP()
# define list of down-regulated probes in new contrast
contrastProbes <- sample(contrast.details[['A2B5 - uns']]$DOWN$ID,50) # random selection of 50 genes
contrast.details <- findDOWN()
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
