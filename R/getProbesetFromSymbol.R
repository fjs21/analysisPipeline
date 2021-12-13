getProbesetFromSymbol <-
function(syms) {
  require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
  ps <- unlist(mget(as.character(syms), revmap(eval(parse(text = paste(chipAnnotation,"SYMBOL",sep="")))), ifnotfound = NA))
  ps <- ps[!is.na(ps)]
  return(ps)
 }

