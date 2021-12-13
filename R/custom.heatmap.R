custom.heatmap <-
function(heatmap.probes,
 eset, title = "Heatmap",
 groups = NULL, maxLabels = 50, sampleLabels = NULL, geneLabels = TRUE,
 relativeTo = "all",
 fnHeatmapCol = heatmapCol, lim = NULL,
 dendrogram = "row",
 Colv = FALSE, ... ) {

 if (is(eset, "ExpressionSet")) {
  #data = exprs(eset)[featureNames(eset) %in% heatmap.probes,]
  probes = featureNames(eset)[match(heatmap.probes, featureNames(eset))]
  probes = unique(probes[!is.na(probes)])
  data = exprs(eset)[probes,]
  if (!exists("chipAnnotation")) chipAnnotation <<- annotation(eset)
  if (is.null(groups)) { groups = eset$phenotype }
 } else {
  if (class(eset) == "matrix") {
   #data = eset[rownames(eset) %in% heatmap.probes,]	
   probes = rownames(eset)[match(heatmap.probes, rownames(eset))]
   probes = unique(probes[!is.na(probes)])
   data = eset[probes,]
   if (is.null(groups)) {
    cat("No groups supplied, using column names\n")
    groups = factor(colnames(eset)) }
   if (!exists("chipAnnotation")) {
    stop("chipAnnotation must be set for data supplied as a matrix.")
   }
  }
  else {
   stop("customTree currently only supports ExpressionSet and matrices")
  }
 } 
 
 suppressWarnings(if (!is.null(relativeTo)) {
  if (relativeTo == "all") { 
   data.medians = apply(data,1,median)
   data <- sweep(data, 1, data.medians)
  } else {
   if (relativeTo == "none") {
   
   } else {
    data.medians = apply(data,1,function(x) median(x[relativeTo]))
    data <- sweep(data, 1, data.medians)
   }
  }
 })

 sample.order <- order(groups) # uses levels from groups
 phenotype.groups <- as.integer(groups[sample.order])
 if(max(phenotype.groups) > 12) { 
  csc <- rainbow(max(phenotype.groups))[phenotype.groups]
 } else {
  csc <- suppressWarnings(brewer.pal(max(phenotype.groups),"Set3")[phenotype.groups])
 }

 hmcol = colorRampPalette(brewer.pal(10, "RdYlGn")[10:1])(256)
 hmcol2 = fnHeatmapCol(data = data, col = hmcol, lim = lim) 

 if(dim(data)[1] > maxLabels) {labRow = ""} else {
  if(geneLabels) {
   labRow = getSymbolFromProbeset(rownames(data))
  } else {
   labRow = rownames(data)
  } 
 }

 if(is.null(sampleLabels)) sampleLabels <- groups

 heatmap.2(data[,sample.order],
  dendrogram = dendrogram,
  Colv = Colv,
  ColSideColors = csc, 
  col = hmcol2,
  labCol = sampleLabels[sample.order],
  labRow = labRow,
  margins = c(8,10),
  trace = "none",
  density.info = "none",
  lhei = c(1,4),
  main = title,
# new for R v3.1.1
  symbreaks = F, symkey = F,
  ... )

 return(data)
}

