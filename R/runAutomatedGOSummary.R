runAutomatedGOSummary <-
function(contrast.details) {
 # code requires that chipAnnotation is set
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

 GOSummaryUP <- function(i) {
  contrast <- names(contrast.details)[i]
  write.csv(GOSummary(contrast.details[[i]]$UPau), file = paste(contrast,"_UP_GOsummary.csv", sep=""))
  cat("Done UP genes.\n")
 }

 GOSummaryDOWN <- function(i) {
  contrast <- names(contrast.details)[i]
  write.csv(GOSummary(contrast.details[[i]]$DOWNau), file = paste(contrast,"_DOWN_GOsummary.csv", sep=""))
  cat("Done DOWN genes\n")
 }

 for (i in 1:length(names(contrast.details))) {
  contrast <- names(contrast.details)[i]
  cat(contrast,"\n")
  tryCatch(GOSummaryUP(i), error = function(e) cat("Error with UP genes in",contrast,".\n"))
  tryCatch(GOSummaryDOWN(i), error = function(e) cat("Error with DOWN genes in",contrast,".\n"))
 }

}

