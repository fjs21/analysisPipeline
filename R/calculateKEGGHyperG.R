calculateKEGGHyperG <-
function(contrast, GeneIDUniverse, contrast.details) {
 require(annaffy)
 require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
 require(GOstats)

 contrast.hyperG <- list() 

 if(length(contrast.details[[contrast]]$UPa$ID) != 0) {
  geneIds = unique(unlist(getGeneIDFromProbeset(contrast.details[[contrast]]$UPa$ID)))
  params = new("KEGGHyperGParams",
   geneIds = as.numeric(geneIds),
   universeGeneIds = as.numeric(GeneIDUniverse), annotation = chipAnnotation,
   pvalueCutoff = 0.05, 
   testDirection = "over")
  tryCatch({
   contrast.hyperG$UP_KEGG <- hyperGTest(params)
   cat(dim(summary(contrast.hyperG$UP_KEGG,p=0.05))[1],"KEGG significant pathways found in UP genes.\n")
   write.csv(summary(contrast.hyperG$UP_KEGG,p=0.05), file=paste(contrast,"UP_KEGGHyperG.csv",sep="_"))
  }, error = function(e) cat("Error occured in HyperG.\n")) 
 }

 if(length(contrast.details[[contrast]]$DOWNa$ID) != 0) {
  geneIds = unique(unlist(getGeneIDFromProbeset(contrast.details[[contrast]]$DOWNa$ID)))
  params = new("KEGGHyperGParams",
   geneIds = as.numeric(geneIds),
   universeGeneIds = as.numeric(GeneIDUniverse), annotation = chipAnnotation,
   pvalueCutoff = 0.05, 
   testDirection = "over")
  tryCatch({
   contrast.hyperG$DOWN_KEGG <- hyperGTest(params)
   cat(dim(summary(contrast.hyperG$DOWN_KEGG,p=0.05))[1],"KEGG significant pathways found in DOWN genes.\n")
   write.csv(summary(contrast.hyperG$DOWN_KEGG, p=0.05), file=paste(contrast,"DOWN_KEGGHyperG.csv",sep="_"))
  }, error = function(e) cat("Error occured in HyperG.\n")) 
 }
 return(contrast.hyperG)
}

