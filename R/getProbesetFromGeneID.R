getProbesetFromGeneID <-
function(EG) {
  require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
  ps <- unlist(mget(as.character(EG), revmap(eval(parse(text = paste(chipAnnotation,"ENTREZID",sep="")))), ifnotfound = NA))
  ps <- ps[!is.na(ps)]
  return(ps)
 }

