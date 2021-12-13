getDescriptionFromProbeset <-
function(probes) {
  require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
  #sapply(probes, function(y) {
  # desc <- as.character(get(y, eval(parse(text = paste(chipAnnotation,"GENENAME",sep="")))))
  # if (length(desc) == 0) "NA" else desc	
  # })
  desc <- mget(probes,eval(parse(text = paste(chipAnnotation,"GENENAME", sep = ""))), ifnotfound = NA)
  desc <- sapply(desc, function(x) as.character(x[1]))
  desc[is.na(desc)] <- "NA"
  return(desc)
}

