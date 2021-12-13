convertHumanToBovineGeneID <-
function(HsEG, verbose = F) {
 # Check for valid input
 if (length(HsEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 HsEG <- as.character(HsEG)

 cat("Converting human geneIds to porcine.")
 export = rep(0, length(HsEG))
 names(export) = HsEG  
 cat(".") 

# Homologene based
 require(annotationTools)
 if (exists("homologene")) {
  suppressWarnings(SsEG <- getHOMOLOG(HsEG, 9913, homologene))
  SsEG <- sapply(SsEG, function(x) {as.character(x[1])}) 
  names(SsEG) = HsEG
  SsEG = SsEG[!is.na(SsEG)]
  export[names(export) %in% names(SsEG)] = SsEG
  if(verbose) cat(sum(names(export) %in% names(SsEG)))
  cat(".")
 } else {
  warning("Homologene is not loaded!\n")
 }

# Ensembl BioMart based
 if (any(export == 0)) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Ss.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   HsENSG = getBM(attributes = c("entrezgene","ensembl_gene_id"), filters = "entrezgene",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   SsENSG = getBM(attributes = c("ensembl_gene_id","TODO_homolog_ensembl_gene"), filters = "entrezgene",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(SsENSG)[1] == 0) stop("No SsENSG.")
   SsEG = getBM(attributes =  c("ensembl_gene_id","entrezgene"), filters = "ensembl_gene_id",
    values = SsENSG[,2], mart = Ss.ensembl)
   SsEG = SsEG[!is.na(SsEG$entrezgene),]
   # match up Human to Mouse results
   SsENSG$SsEG = SsEG[match(SsENSG$btaurus_homolog_ensembl_gene, SsEG$ensembl_gene_id),2]
   HsENSG$SsEG = SsENSG[match(HsENSG$ensembl_gene_id, SsENSG$ensembl_gene_id),3]
   SsEG = as.character(HsENSG$SsEG)
   names(SsEG) = HsENSG[,1]
   SsEG = SsEG[!is.na(SsEG)]
   SsEG = SsEG[unique(names(SsEG))]
   # export data
   for (i in which(names(export) %in% names(SsEG))) {
    export[i] <- SsEG[names(SsEG) == names(export[i])] 
   }
   if(verbose) cat(length(SsEG))
   cat(".")
  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 }

# Bioconductor-based
 if (any(export == 0)) {
  require("org.Hs.eg.db")
  HsENSP = mget(HsEG[export == 0], org.Hs.egENSEMBLPROT, ifnotfound = NA)
  HsENSP = sapply(HsENSP[!is.na(HsENSP)], function(x) {x[1]}) # select first match only 
  require("hom.Hs.inp.db")
  tryCatch({
   if (length(HsENSP) == 0) stop("No HsENSP.")
   SsENSP = mget(HsENSP, hom.Hs.inpBOSTA, ifnotfound = NA)
   SsENSP = sapply(SsENSP[!is.na(SsENSP)], function(x) {x[1]}) # select first match only
   if (length(SsENSP) == 0) stop("No SsENSP.")   
   require("org.Ss.eg.db")
   SsEG = mget(SsENSP, org.Ss.egENSEMBLPROT2EG, ifnotfound = NA)
   SsEG = sapply(SsEG, function(x) {x[1]}) # select first match only
   SsEG = SsEG[!is.na(SsEG)]
   # match up Human to Bovine results 
   names(SsEG) = names(HsENSP)[HsENSP %in% names(SsENSP)[SsENSP %in% names(SsEG)]]
   SsEG = SsEG[!is.na(names(SsEG))]
   # prepare data export
   export[names(export) %in% names(SsEG)] = SsEG
   if(verbose) cat(sum(names(export) %in% names(SsEG)))
   cat(".")
  }, error = function(e) ifelse(verbose,print(e),cat("e")))
 }
 
 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 return(export)
}

