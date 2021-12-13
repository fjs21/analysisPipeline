marker.heatmap <-
function(eset, 
 marker.probes = get("marker.probes", envir = environment(marker.heatmap)),
 clustered = FALSE, title = "Heatmap - Marker Genes",
 groups = NULL, maxLabels = 50, sampleLabels = NULL, geneLabels = TRUE,
 relativeTo = "all", fnHeatmapCol = heatmapCol, lim = NULL ) {

 if (exists('markerGenes')) {
  if (length(markerGenes$SYMBOL) == 0)
   stop("markerGenes must contain a 'SYMBOL' column.")
  if (length(markerGenes$TblOrder) == 0)
   stop("markerGenes must contain a 'TblOrder' column.")
 } else {
  stop("markerGenes object does not exist, please define.")
 }

 if (is(eset, "ExpressionSet")) {
  data = exprs(eset)[featureNames(eset) %in% marker.probes,]
  if (is.null(groups)) { groups = eset$phenotype }
 } else {
  if (class(eset) == "matrix") {
   data = eset[rownames(eset) %in% marker.probes,]     
   if (is.null(groups)) {
    cat("No groups supplied, using column names\n")
    groups = factor(colnames(eset)) }
  }
  else {
   stop("marker.heatmap currently only supports ExpressionSet and matrices")
  }
 } 

 if (relativeTo != 0) {
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
 }

 if (dim(data)[1] == 0) stop("no matching marker.probes in eset variable")

 sample.order <- order(groups) # uses levels from groups
 phenotype.groups <- as.integer(groups[sample.order])
 if(max(phenotype.groups) > 12) { 
  csc <- rainbow(max(phenotype.groups))[phenotype.groups]
 } else {
  csc <- suppressWarnings(brewer.pal(max(phenotype.groups),"Set3")[phenotype.groups])
 }

 hmcol = colorRampPalette(brewer.pal(10, "RdYlGn")[10:1])(256)
 hmcol2 = fnHeatmapCol(data = data, col = hmcol, lim = lim) 

 if(geneLabels) {
  labRow = getSymbolFromProbeset(rownames(data))
 } else {
  labRow = rownames(data)
 } 
 
 if(is.null(sampleLabels)) sampleLabels <- groups

 # add colors by cell type (find GeneIDs then CellType info)
 TblPosition <- factor(sapply(labRow, 
  function(x) markerGenes$TblOrder[markerGenes$SYMBOL %in% x][1]))
 gene.order <- order(TblPosition, decreasing = FALSE)
 cellTypes <- as.integer(TblPosition[gene.order]) 
 rsc <- suppressWarnings(brewer.pal(max(cellTypes), "Accent")[cellTypes])

 suppressWarnings(heatmap.2(data[gene.order,sample.order], 
  dendrogram = "both",
  Rowv = clustered,
  Colv = clustered,
  RowSideColors = rsc,
  ColSideColors = csc, 
  col = hmcol2,
  labCol = sampleLabels[sample.order],
  labRow = labRow[gene.order],
  margins = c(8,10),
  trace = "none",
  density.info = "none",
  main = title,
  # new for R v3.1.1
  symbreaks = F, symkey = F))
}

