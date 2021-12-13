parseContrasts <-
function(fit3, FCcutoff = 2, FDRcutoff = 0.05, adjust = "fdr") {
 require(limma)
 contrast.details <- list()
 for (i in 1:length(colnames(fit3$contrasts))) {
	cat(colnames(fit3$contrasts)[i],"\n")
	contrast <- colnames(fit3$contrasts)[i]
	contrast.details[[contrast]]$topT <- 
		topTable(fit3, coef=i, p.value=FDRcutoff, lfc=log2(FCcutoff), number=100000, adjust.method = adjust)
	contrast.details[[contrast]]$topT$ID <- rownames(contrast.details[[contrast]]$topT)

     # Find all probes UP
	contrast.details[[contrast]]$UP <-
		contrast.details[[contrast]]$topT[contrast.details[[contrast]]$topT$logFC > 0,]
	if(dim(contrast.details[[contrast]]$UP)[1] == 0) {cat("No significant probe sets UP regulated.\n")} else {
	 cat("Probe sets UP -",as.character(length(contrast.details[[contrast]]$UP[,1])),"\n")
	# Annotate probes UP	
	contrast.details[[contrast]]$UPa <- 
		data.frame(Symbol = getSymbolFromProbeset(contrast.details[[contrast]]$UP$ID), 
			Description = getDescriptionFromProbeset(contrast.details[[contrast]]$UP$ID),
			GeneID = getGeneIDFromProbeset(contrast.details[[contrast]]$UP$ID),
			contrast.details[[contrast]]$UP)
	write.csv(contrast.details[[contrast]]$UPa, 
		file=paste(contrast, "_UP.csv", sep=""))
	# Find unique genes UP
	contrast.details[[contrast]]$UPau <- 
		contrast.details[[contrast]]$UPa[!duplicated(contrast.details[[contrast]]$UPa$GeneID),]
	contrast.details[[contrast]]$UPau <- contrast.details[[contrast]]$UPau[contrast.details[[contrast]]$UPau$GeneID != 0,]
	cat("Unique genes UP -",as.character(dim(contrast.details[[contrast]]$UPau)[1]),"\n")
	write.csv(contrast.details[[contrast]]$UPau, 
		file=paste(contrast, "_UP_uniques.csv", sep=""))
	if(exists('markerGenes')) { write.csv(contrast.details[[contrast]]$UPau[contrast.details[[contrast]]$UPau$GeneID %in% markerGenes$GeneID,],file=paste(contrast, "_UP_markers.csv", sep="")) }
	}

	# Find all probes DOWN
	contrast.details[[contrast]]$DOWN <-
		contrast.details[[contrast]]$topT[contrast.details[[contrast]]$topT$logFC < 0,]
	if(dim(contrast.details[[contrast]]$DOWN)[1] == 0) {cat("No significant probe sets DOWN regulated.\n")} else {
	cat("Probe sets DOWN -",as.character(length(contrast.details[[contrast]]$DOWN[,1])),"\n")
	# Annotate probes DOWN
	contrast.details[[contrast]]$DOWNa <- 
		data.frame(Symbol = getSymbolFromProbeset(contrast.details[[contrast]]$DOWN$ID), 
			Description = getDescriptionFromProbeset(contrast.details[[contrast]]$DOWN$ID),
			GeneID = getGeneIDFromProbeset(contrast.details[[contrast]]$DOWN$ID),
			contrast.details[[contrast]]$DOWN)
	write.csv(contrast.details[[contrast]]$DOWNa, 
		file=paste(contrast, "_DOWN.csv", sep=""))
	# Find unique genes DOWN
	contrast.details[[contrast]]$DOWNau <- 
		contrast.details[[contrast]]$DOWNa[!duplicated(contrast.details[[contrast]]$DOWNa$GeneID),]
	contrast.details[[contrast]]$DOWNau <- contrast.details[[contrast]]$DOWNau[contrast.details[[contrast]]$DOWNau$GeneID != 0,]
	cat("Unique genes DOWN -",as.character(dim(contrast.details[[contrast]]$DOWNau)[1]),"\n")
	write.csv(contrast.details[[contrast]]$DOWNau, 
		file=paste(contrast, "_DOWN_uniques.csv", sep=""))
	if(exists('markerGenes')) {write.csv(contrast.details[[contrast]]$DOWNau[contrast.details[[contrast]]$DOWNau$GeneID %in% markerGenes$GeneID,],file=paste(contrast, "_DOWN_markers.csv", sep=""))}
	}
 }
 return(contrast.details)
}

