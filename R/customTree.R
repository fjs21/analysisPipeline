customTree <-
function(eset, hc_labels = NULL, groups = NULL, title = "Hierarchical Clustering" ) {

 if (is(eset, "ExpressionSet")) {
  e = exprs(eset)
  if(is.null(hc_labels)) hc_labels = sampleNames(eset)
  if(is.null(groups) & !is.null(eset$phenotype)) groups = eset$phenotype
 } else {
  if (class(eset) == "matrix") {
   e = eset
   if(is.null(hc_labels)) hc_labels = colnames(eset)
   if(is.null(groups)) groups <- factor(hc_labels) 
  }
  else {
   stop("customTree currently only supports ExpressionSet and matrices")
  }
 }
 hc = hclust(dist(t(e), method="euclidean")) 
 hc$labels = paste(hc_labels,1:length(hc_labels),sep="__")
 dend = as.dendrogram(hc, hang = 0.1)
 if(length(levels(groups)) > 8) { 
  treeCol <- rainbow(length(levels(groups)))
 } else {
  treeCol <- suppressWarnings(brewer.pal(length(levels(groups)),"Dark2"))
 } 
 dend = dendrapply(dend, function(n) {
   a <- attributes(n)
   if(is.leaf(n)) {
     attr(n, "nodePar") <- list(pch = NA, lab.col = treeCol[as.numeric(groups[which(hc$labels == a$label)])])
     attr(n, "edgePar") <- list(col = treeCol[as.numeric(groups[which(hc$labels == a$label)])])
     attr(n, "label") <- substring(attr(n, "label"),1,regexpr("__",attr(n, "label"))-1)
   } else {
    nc = dendrapply(n, function(subn) {
     a <- attributes(subn)
     as.character(a$label)
    })
    nc = unlist(nc)
    nc = nc[!is.na(nc)]
    nc = groups[hc$labels %in% nc]
    if(length(unique(nc)) == 1) {
     attr(n, "edgePar") <- list(col = treeCol[as.numeric(nc)])    
    } else {
     attr(n, "edgePar") <- list(col = 1)
    }
   }
   n
  })
  plot(dend, main = title)
  legend("bottom", legend = levels(groups),lwd = 1, col = treeCol, ncol = 2)
 }

