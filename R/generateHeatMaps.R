generateHeatMaps <-
function(eset, export = TRUE, contrast.names = NULL, ...) {
 if(is.null(contrast.names))  contrast.names <- names(contrast.details)
 for (coef in 1:length(contrast.names)) {
  if(export) { pdf(file = paste( contrast.names[coef],"UP.pdf", sep="_")) 
  } else { par(ask=TRUE) }
  tryCatch(
   custom.heatmap(contrast.details[[coef]]$UPau$ID[1:50], eset,
    title = paste(contrast.names[coef],"UP"), ...)
  ,error = function(e) cat("Could not generate heatmap of genes UP in contrast =", contrast.names[coef],"\n") )
  if(export) { dev.off()
   pdf(file = paste(contrast.names[coef],"DOWN.pdf", sep="_"))}
  tryCatch(
   custom.heatmap(contrast.details[[coef]]$DOWNau$ID[1:50], eset,
    title = paste(contrast.names[coef],"DOWN"), ...)
  ,error = function(e) cat("Could not generate heatmap of genes DOWN in contrast =", contrast.names[coef],"\n") )
  if(export) { dev.off() } else { par(ask=FALSE) }
 }
}

