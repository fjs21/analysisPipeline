findDOWN <-
function() {
 cat(contrast,"\n")
 contrast.details[[contrast]]$DOWN = contrast.details[[contrastRatio]]$DOWN[contrast.details[[contrastRatio]]$DOWN$ID %in% 
  contrastProbes,]
 cat("Probe sets DOWN -",as.character(length(contrast.details[[contrast]]$DOWN[,1])),"\n")
 contrast.details[[contrast]]$DOWNa <- 
  data.frame(Symbol = getSymbolFromProbeset(contrast.details[[contrast]]$DOWN$ID),
  Description = getDescriptionFromProbeset(contrast.details[[contrast]]$DOWN$ID), 
  GeneID = getGeneIDFromProbeset(contrast.details[[contrast]]$DOWN$ID),
  contrast.details[[contrast]]$DOWN)
 # Write file
 write.csv(contrast.details[[contrast]]$DOWNa, 
  file=paste(contrast, "_DOWN.csv", sep=""))
 # Find unique genes DOWN
 contrast.details[[contrast]]$DOWNau <- 
  contrast.details[[contrast]]$DOWNa[!duplicated(contrast.details[[contrast]]$DOWNa$GeneID),]
 contrast.details[[contrast]]$DOWNau <- contrast.details[[contrast]]$DOWNau[contrast.details[[contrast]]$DOWNau$GeneID != 0,]
 cat("Unique genes DOWN -",as.character(dim(contrast.details[[contrast]]$DOWNau)[1]),"\n")
 write.csv(contrast.details[[contrast]]$DOWNau, 
  file=paste(contrast, "_DOWN_uniques.csv", sep=""))
 if(exists('markerGenes')) {write.csv(contrast.details[[contrast]]$DOWNau[contrast.details[[contrast]]$DOWNau$GeneID %in% markerGenes$GeneID,],
  file=paste(contrast, "_DOWN_markers.csv", sep=""))}
 return(contrast.details)
 }

