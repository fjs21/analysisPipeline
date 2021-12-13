getSymbolFromProbeset <-
function(probes) {	
  require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
#  sapply(probes, function(y) {
#   sym <- as.character(get(y, eval(parse(text = paste(chipAnnotation,"SYMBOL",sep="")))))
#   if (length(sym) == 0) "NA" else sym	
#   })
  syms <- mget(probes,eval(parse(text = paste(chipAnnotation,"SYMBOL", sep = ""))), ifnotfound = NA)
  syms <- sapply(syms, function(x) as.character(x[1]))
  syms[is.na(syms)] <- 0
  return(syms)
  }

