plot3dPCA <-
function (eset, groups = NULL, groupnames = NULL, addtext = NULL, 
    x.coord = NULL, y.coord = NULL, screeplot = FALSE, squarepca = FALSE, 
    pch = NULL, col = NULL, pcs = c(1, 2, 3), legend = TRUE, ...) 
{
 require(lattice)
 require(rgl)

 if (length(pcs) != 3)
  stop("You can only plot three principal components.\n", call. = FALSE)
  # ExpressionSet object
  if (is(eset, "ExpressionSet")) {
   if (max(pcs) > dim(exprs(eset))[2]) 
    stop(paste("There are only", dim(exprs(eset))[2], 
     "principal components to plot.\n", call. = FALSE))
   if (is.null(groups))       groups = factor(pData(eset)$phenotype)
   if (is.null(groupnames))   groupnames <- levels(groups)
   if (is.factor(groupnames)) groupnames <- as.character(groupnames)
   pca <- prcomp(t(exprs(eset)))
   len <- length(sampleNames(eset))
   sampleNames <- sampleNames(eset)
  } else {
   # Matrix object
   if (class(eset) == "matrix") {
    if (max(pcs) > dim(eset)[2]) 
     stop(paste("There are only", dim(eset)[2], "principal components to plot.\n", 
      call. = FALSE))
    if (is.null(groupnames) & is.null(groups)) {
     groupnames <- colnames(eset) 
    } else {
     if (!is.null(groups)) groupnames <- levels(groups)
    }
    if (is.factor(groupnames)) groupnames <- as.character(groupnames)
    pca <- prcomp(t(eset))
    len <- dim(eset)[2]
    sampleNames = colnames(eset)
   } else {
    stop("plotPCA currently only supports ExpressionSet and matrices")
   }
  }
  if (screeplot) {
   plot(pca, main = "Screeplot")
  } else {   
   x = pca$x[,pcs[1]]
   y = pca$x[,pcs[2]]
   z = pca$x[,pcs[3]]
   if (squarepca) {
    ylim <- max(abs(range(x)))
    ylim <- c(-ylim, ylim)
    zlim <- ylim
   } else {
    ylim <- range(y, finite = TRUE)
    zlim <- range(z, finite = TRUE)
   }
   # With groups defined - color and shape of spots based on groups
   if (!is.null(groups)) {
    if (is.null(pch)) pch <- as.numeric(groups) - 1
    if (is.null(col)) {
     if (length(levels(groups)) > 8) {
      col <- rainbow(length(levels(groups)))[as.numeric(groups)]
     } else {
      col <- suppressWarnings(brewer.pal(length(levels(groups)),"Dark2")[as.numeric(groups)])
     }	
    }
    plot3d(x,y,z, 
     zlab = paste("PC", pcs[3], sep = ""),
     ylab = paste("PC", pcs[2], sep = ""), 
     xlab = paste("PC", pcs[1], sep = ""), 
     main = "Principal Components Plot", 
     zlim = zlim, ylim = ylim, col = col, 
     type = "s", radius = 2,  box = FALSE, 
     ...)
    if (is.null(addtext)) {
   } else {
    text3d(x, y+5, z+5, addtext, col = col, cex = 0.7)
   }
   if (legend) {
    print(cloud(z ~ x * y, 
     zlab = paste("PC", pcs[3], sep = ""), 
     ylab = paste("PC", pcs[2], sep = ""), 
     xlab = paste("PC", pcs[1], sep = ""), 
     main = "Principal Components Plot", 
     zlim = zlim, ylim = ylim, pch = pch, col = col,
     R.mat = par3d("userMatrix"), screen = list(),
     key = list(corner = c(0,1), border = TRUE, 
      points = list(pch=0:(length(levels(groups)) - 1),
       col = if (length(levels(groups)) > 8) {
        rainbow(length(levels(groups)))
       } else {
        suppressWarnings(brewer.pal(8,"Dark2"))[1:length(levels(groups))]
       }),
      text = list(groupnames))
     ,...)) 
    } else {
    print(cloud(z ~ x * y, 
     zlab = paste("PC", pcs[3], sep = ""), 
     ylab = paste("PC", pcs[2], sep = ""), 
     xlab = paste("PC", pcs[1], sep = ""), 
     main = "Principal Components Plot", 
     zlim = zlim, ylim = ylim, pch = pch, col = col,
     R.mat = par3d("userMatrix"), screen = list(),
     ,...)) 
   } 
  } else {
   # No groups defined
   if (is.null(pch)) pch <- 0:(len - 1)
   if (is.null(col)) col <- 1:len
   print(plot3d(x,y,z, 
    zlab = paste("PC", pcs[3], sep = ""), 
    ylab = paste("PC", pcs[2], sep = ""), 
    xlab = paste("PC", pcs[1], sep = ""),
    main = "Principal Components Plot", 
    zlim = zlim, ylim = ylim, col = col, 
    type = "s", radius = 2, box = FALSE,
    ...))
   if (is.null(addtext)) {
   } else {
    text3d(x, y+5, z+5, addtext, col = col, cex = 0.7)
   }
   if (legend) {
    print(cloud(z ~ x * y, 
     zlab = paste("PC", pcs[3], sep = ""), 
     ylab = paste("PC", pcs[2], sep = ""), 
     xlab = paste("PC", pcs[1], sep = ""), 
     main = "Principal Components Plot", 
     zlim = zlim, ylim = ylim, pch = pch, col = col,
     R.mat = par3d("userMatrix"), screen = list(),
     key = list(corner = c(0,1), border = TRUE, 
      points = list(pch=pch, col=col),
      text = list(sampleNames))
     ,...))
   } else {
    print(cloud(z ~ x * y, 
     zlab = paste("PC", pcs[3], sep = ""), 
     ylab = paste("PC", pcs[2], sep = ""), 
     xlab = paste("PC", pcs[1], sep = ""), 
     main = "Principal Components Plot", 
     zlim = zlim, ylim = ylim, pch = pch, col = col,
     R.mat = par3d("userMatrix"), screen = list()
     ,...))
   }
  }
 }
}