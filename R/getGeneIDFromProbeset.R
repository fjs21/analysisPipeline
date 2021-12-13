getGeneIDFromProbeset <-
function(probes) {
  require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
  #sapply(probes, function(y) {
  # GeneID = as.numeric(get(y, eval(parse(text = paste(chipAnnotation,"ENTREZID",sep="")))))
  # GeneID[is.na(GeneID)] = 0
  # if (length(GeneID) == 0) 0 else GeneID
  #})

  EG <- mget(probes,eval(parse(text = paste(chipAnnotation,"ENTREZID", sep = ""))), ifnotfound = NA)
  EG <- sapply(EG, function(x) as.numeric(x[1]))
  EG[is.na(EG)] <- 0
  return(EG)
}

