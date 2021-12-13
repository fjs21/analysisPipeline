findUP <-
function() {
 cat(contrast,"\n")
 contrast.details[[contrast]]$UP = contrast.details[[contrastRatio]]$UP[contrast.details[[contrastRatio]]$UP$ID %in% 
  contrastProbes,]
 cat("Probe sets UP -",as.character(length(contrast.details[[contrast]]$UP[,1])),"\n")
 contrast.details[[contrast]]$UPa <- 
  data.frame(Symbol = getSymbolFromProbeset(contrast.details[[contrast]]$UP$ID),
  Description = getDescriptionFromProbeset(contrast.details[[contrast]]$UP$ID), 
  GeneID = getGeneIDFromProbeset(contrast.details[[contrast]]$UP$ID),
  contrast.details[[contrast]]$UP)
 # Write file
 write.csv(contrast.details[[contrast]]$UPa, 
  file=paste(contrast, "_UP.csv", sep=""))
 # Find unique genes UP
 contrast.details[[contrast]]$UPau <- 
  contrast.details[[contrast]]$UPa[!duplicated(contrast.details[[contrast]]$UPa$GeneID),]
 contrast.details[[contrast]]$UPau <- contrast.details[[contrast]]$UPau[contrast.details[[contrast]]$UPau$GeneID != 0,]
 cat("Unique genes UP -",as.character(dim(contrast.details[[contrast]]$UPau)[1]),"\n")
 write.csv(contrast.details[[contrast]]$UPau, 
  file=paste(contrast, "_UP_uniques.csv", sep=""))
 if(exists('markerGenes')) {write.csv(contrast.details[[contrast]]$UPau[contrast.details[[contrast]]$UPau$GeneID %in% markerGenes$GeneID,],
  file=paste(contrast, "_UP_markers.csv", sep=""))}
 return(contrast.details)
 }

