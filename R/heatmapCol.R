heatmapCol <- function (data, col, lim = NULL) {
    nrcol <- length(col)
    data.range <- range(data)
    if (diff(data.range) == 0)
        stop("data has range 0")
    if(is.null(lim)) 
        lim <- min(abs(data.range))
    if (lim <= 0)
        stop("lim has to be positive")
    if (lim > min(abs(data.range))) {
        warning("specified bound 'lim' is out of data range\n\nhence 'min(abs(range(data)))' is used")
        lim <- min(abs(data.range))
    }
    nrcol <- length(col)
    reps1 <- ceiling(nrcol * (-lim - data.range[1])/(2 * lim))
    reps2 <- ceiling(nrcol * (data.range[2] - lim)/(2 * lim))
    col1 <- c(rep(col[1], reps1), col, rep(col[nrcol], reps2))
    return(col1)
}