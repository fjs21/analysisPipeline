calculateGOHyperG <-
function(contrast, eset, contrast.details) {
 contrast.hyperG <- list()
 if(length(contrast.details[[contrast]]$UPa$ID) != 0) {
  contrast.hyperG$UP_GO <- analysisPipeline::topGO(contrast.details[[contrast]]$UPa$ID, eset)
  write.csv(contrast.hyperG$UP_GO, file=paste(contrast,"UP_topGO.csv",sep="_"))
 }
 if(length(contrast.details[[contrast]]$DOWNa$ID) != 0) {
  contrast.hyperG$DOWN_GO <- analysisPipeline::topGO(contrast.details[[contrast]]$DOWNa$ID, eset)
  write.csv(contrast.hyperG$DOWN_GO, file=paste(contrast,"DOWN_topGO.csv",sep="_"))
 }
 return(contrast.hyperG)
}

