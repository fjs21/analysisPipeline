heatmapCol_stretch <- function(data, col, lim = NULL) {
 nrcol <- length(col)
 data.range <- range(data)
 if (diff(data.range) == 0)
        stop("data has range 0")
 if(is.null(lim)) {
  lim <- max(abs(data.range))
 } else {
  if(lim>max(abs(data.range))) {
   warning("lim value too high")
   lim <- max(abs(data.range))
  }
 }

 nrup <- data.range[2]
 nrdown <- -data.range[1]
 nr <- nrup + nrdown
 nrdown <- ceiling(nrdown/nr * nrcol)
 nrup <- ceiling(nrup/nr * nrcol)

 midP <- ceiling(nrcol/2)

 if(lim < abs(data.range[1])) {
   col1 <- rep(1, ceiling(nrdown - (lim/-data.range[1] * nrdown)))
   col2 <- ceiling(seq(1,midP, length.out = lim/-data.range[1] * nrdown))
 } else {
   col1 <- NULL
   col2 <- ceiling(seq(1,midP, length.out = nrdown))
 }
 if(lim < abs(data.range[2])) { 
   col3 <- ceiling(seq(midP,nrcol, length.out = lim/data.range[2] * nrup))
   col4 <- rep(nrcol,ceiling(nrup - (lim/data.range[2] * nrup)))
 } else {
   col3 <- ceiling(seq(midP,nrcol, length.out = nrup))
   col4 <- NULL
 }
 return(col[c(col1,col2[2:length(col2)-1],col3[1:length(col3)-1],col4)])
}
