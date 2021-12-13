topGO <-
function(probes, eset = eset) {

 if (!exists('mapping')) {mapping = 'org.Hs.eg'}
 if (!exists('chipAnnotation')) {chipAnnotation = annotation(eset)}
 
 require(topGO)
 require(paste(mapping,".db",sep=""),character.only = TRUE)
 require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)

 geneIds = unique(unlist(mget(probes, 
  eval(parse(text = paste(chipAnnotation,"ENTREZID",sep=""))), ifnotfound = NA )))

 allGeneNames = unique(unlist(mget(featureNames(eset), 
  eval(parse(text = paste(chipAnnotation,"ENTREZID",sep=""))), ifnotfound = NA )))
 allGenes <- factor(as.integer(allGeneNames %in% geneIds))
 names(allGenes) = allGeneNames

 GOdata <- new("topGOdata", ontology = "BP", allGenes = allGenes,
  annotationFun = annFUN.org, mapping = mapping)

 elim.fisher <- runTest(GOdata, algorithm = "elim", statistic = "fisher")

 allRes <- GenTable(GOdata, p.value = elim.fisher, orderBy = "p.value", ranksOf = "p.value", topNodes = 100)
 allRes <- allRes[as.numeric(allRes$p.value) < 0.01,]
 allRes$geneIds <- sapply(genesInTerm(GOdata, whichGO = allRes$GO.ID), function(x) {
  x <- x[x %in% geneIds]
  paste(x, collapse="|")
  })
 allRes$symbols <- sapply(genesInTerm(GOdata, whichGO = allRes$GO.ID), function(x) {
  x <- x[x %in% geneIds]
  x <- mget(x, eval(parse(text = paste(mapping,"SYMBOL",sep=""))), ifnotfound = NA )
  paste(x, collapse="|")
  })
 return(allRes)
}

