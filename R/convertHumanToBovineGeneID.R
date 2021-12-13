convertHumanToBovineGeneID <-
function(HsEG, verbose = F) {
 # Check for valid input
 if (length(HsEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 HsEG <- as.character(HsEG)

 cat("Converting human geneIds to bovine.")
 export = rep(0, length(HsEG))
 names(export) = HsEG  
 cat(".") 

# Homologene based
 require(annotationTools)
 if (exists("homologene")) {
  suppressWarnings(BtEG <- getHOMOLOG(HsEG, 9913, homologene))
  BtEG <- sapply(BtEG, function(x) {as.character(x[1])}) 
  names(BtEG) = HsEG
  BtEG = BtEG[!is.na(BtEG)]
  export[names(export) %in% names(BtEG)] = BtEG
  if(verbose) cat(sum(names(export) %in% names(BtEG)))
  cat(".")
 } else {
  warning("Homologene is not loaded!\n")
 }

# Ensembl BioMart based
 if (any(export == 0)) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Bt.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   HsENSG = getBM(attributes = c("entrezgene_id","ensembl_gene_id"), filters = "entrezgene_id",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   BtENSG = getBM(attributes = c("ensembl_gene_id","btaurus_homolog_ensembl_gene"), filters = "entrezgene_id",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(BtENSG)[1] == 0) stop("No BtENSG.")
   BtEG = getBM(attributes =  c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
    values = BtENSG[,2], mart = Bt.ensembl)
   BtEG = BtEG[!is.na(BtEG$entrezgene),]
   # match up Human to Mouse results
   BtENSG$BtEG = BtEG[match(BtENSG$btaurus_homolog_ensembl_gene, BtEG$ensembl_gene_id),2]
   HsENSG$BtEG = BtENSG[match(HsENSG$ensembl_gene_id, BtENSG$ensembl_gene_id),3]
   BtEG = as.character(HsENSG$BtEG)
   names(BtEG) = HsENSG[,1]
   BtEG = BtEG[!is.na(BtEG)]
   BtEG = BtEG[unique(names(BtEG))]
   # export data
   for (i in which(names(export) %in% names(BtEG))) {
    export[i] <- BtEG[names(BtEG) == names(export[i])] 
   }
   if(verbose) cat(length(BtEG))
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
   BtENSP = mget(HsENSP, hom.Hs.inpBOSTA, ifnotfound = NA)
   BtENSP = sapply(BtENSP[!is.na(BtENSP)], function(x) {x[1]}) # select first match only
   if (length(BtENSP) == 0) stop("No BtENSP.")   
   require("org.Bt.eg.db")
   BtEG = mget(BtENSP, org.Bt.egENSEMBLPROT2EG, ifnotfound = NA)
   BtEG = sapply(BtEG, function(x) {x[1]}) # select first match only
   BtEG = BtEG[!is.na(BtEG)]
   # match up Human to Bovine results 
   names(BtEG) = names(HsENSP)[HsENSP %in% names(BtENSP)[BtENSP %in% names(BtEG)]]
   BtEG = BtEG[!is.na(names(BtEG))]
   # prepare data export
   export[names(export) %in% names(BtEG)] = BtEG
   if(verbose) cat(sum(names(export) %in% names(BtEG)))
   cat(".")
  }, error = function(e) ifelse(verbose,print(e),cat("e")))
 }
 
 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 return(export)
}

