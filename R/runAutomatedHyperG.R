runAutomatedHyperG <-
function(eset, contrast.details, saveData = FALSE) {

 # load packages if necessary
 require(annaffy)
 require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
 require(GOstats)

 cat("Starting HyperG analysis...\n")

 GeneIDUniverse = unique(unlist(mget(featureNames(eset), eval(parse(text = paste (chipAnnotation,"ENTREZID",sep=""))), ifnotfound = NA )))
 GeneIDUniverse = GeneIDUniverse[!is.na(GeneIDUniverse)]

 contrast.hyperG <- list()

 # GO based HyperG analysis
 cat("GO Biologic Process HyperG.\n")
 for (i in 1:length(names(contrast.details))) {
  cat(names(contrast.details)[i],"\n")
  contrast <- names(contrast.details)[i]
  contrast.hyperG[[contrast]] <- c(contrast.hyperG[[contrast]],
   calculateGOHyperG(contrast, eset, contrast.details))
 }

 # KEGG based HyperG analysis
 cat("KEGG Pathway HyperG.\n")
 for (i in 1:length(names(contrast.details))) {
  cat(names(contrast.details)[i],"\n")
  contrast <- names(contrast.details)[i]
  contrast.hyperG[[contrast]] <- c(contrast.hyperG[[contrast]],
   calculateKEGGHyperG(contrast, GeneIDUniverse, contrast.details))
 }

 if(saveData) {save.image("HyperG.RData")}

 return(contrast.hyperG)
}

