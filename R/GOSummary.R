GOSummary <-
function(geneTable) {	

 if(!exists('ligands')) {
  require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
  # Ligands
  ligands <<- as.character(get("GO:0005102",
   eval(parse(text = paste(chipAnnotation,"GO2ALLPROBES",sep="")))))
  # Receptors 
  receptors <<- as.character(get("GO:0004872",
   eval(parse(text = paste(chipAnnotation,"GO2ALLPROBES",sep="")))))
  # Transcription Factors 
  TFs <<- as.character(get(c("GO:0006355","GO:0030528", "GO:0003677"),
   eval(parse(text = paste(chipAnnotation,"GO2ALLPROBES",sep="")))))
  # ECM 
  ECM <<- as.character(get(c("GO:0005578", "GO:0007155"),
   eval(parse(text = paste(chipAnnotation,"GO2ALLPROBES",sep="")))))
  # Enzymes
  Enzymes <<- as.character(get(c("GO:0003824","GO:0006096","GO:0016125"),
   eval(parse(text = paste(chipAnnotation,"GO2ALLPROBES",sep="")))))
 }

 cat(dim(geneTable)[1],"genes total.\n")

 if(any(geneTable$ID %in% ligands)) {
  l <- data.frame("Ligands",geneTable[geneTable$ID %in% ligands,])
  names(l) <- c("Type",names(geneTable))
  cat(dim(l)[1],"ligands.\n")
 } else {l = data.frame()}

 if(any(geneTable$ID %in% receptors)) {
  r <- data.frame("Receptors",geneTable[geneTable$ID %in% receptors,])
  names(r) <- c("Type",names(geneTable))
  cat(dim(r)[1],"receptors.\n")
 } else {r = data.frame()}
 
 if(any(geneTable$ID %in% TFs)) {
  tf = data.frame("Transcription Factors",geneTable[geneTable$ID %in% TFs,])
  names(tf) <- c("Type",names(geneTable))
  cat(dim(tf)[1],"transcription factors.\n")
 } else {tf = data.frame()}

 if(any(geneTable$ID %in% ECM)) {
  ecm = data.frame("ECM",geneTable[geneTable$ID %in% ECM,])
  names(ecm) <- c("Type",names(geneTable))
  cat(dim(ecm)[1],"ECM molecules.\n")
 } else {ecm = data.frame()}

 if(any(geneTable$ID %in% Enzymes)) {
  enzymes = data.frame("Enzymes",geneTable[geneTable$ID %in% Enzymes,])
  names(enzymes) <- c("Type",names(geneTable))
  cat(dim(enzymes)[1],"enzymes.\n")
 } else {enzymes = data.frame()}

 Other = !geneTable$ID %in% c(ligands,
  receptors,TFs,ECM,Enzymes)
 other = data.frame("Other",geneTable[Other,])
 names(other) = c("Type",names(geneTable))
 cat(dim(other)[1],"unannotated genes.\n")

 return(rbind(l,r,tf,ecm,enzymes,other))
}

