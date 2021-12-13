convertHumanToRatGeneID <-
function(HsEG, verbose = F) {
 # Check for valid input
 if (length(HsEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 HsEG <- as.character(HsEG)

 cat("Converting human geneIds to rat.")
 export = rep(0, length(HsEG))
 names(export) = HsEG  
 cat(".") 

 # Homologene based (Kuhn et al., 2008 - BMC)
 if (exists("homologene")) { 
  require(annotationTools)
  suppressWarnings(RnEG <- getHOMOLOG(HsEG, 10116, homologene))
  RnEG <- sapply(RnEG, function(x) {as.character(x[1])}) 
  names(RnEG) = HsEG
  RnEG = RnEG[!is.na(RnEG)]
  export[names(export) %in% names(RnEG)] = RnEG
  if(verbose) cat(sum(names(export) %in% names(RnEG)))
  cat(".")
 } else {
  warning("Homologene is not loaded!\n")
 }

 # Ensembl BioMart based
 if (any(export == 0)) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Rn.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   HsENSG = getBM(attributes = c("entrezgene_id","ensembl_gene_id"), filters = "entrezgene_id",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   RnENSG = getBM(attributes = c("ensembl_gene_id","rnorvegicus_homolog_ensembl_gene"), filters = "entrezgene_id",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(RnENSG)[1] == 0) stop("No RnENSG.")
   RnEG = getBM(attributes =  c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
    values = RnENSG[,2], mart = Rn.ensembl)
   RnEG = RnEG[!is.na(RnEG$entrezgene),]
   # match up Human to Mouse results
   RnENSG$RnEG = RnEG[match(RnENSG$rnorvegicus_homolog_ensembl_gene, RnEG$ensembl_gene_id),2]
   HsENSG$RnEG = RnENSG[match(HsENSG$ensembl_gene_id, RnENSG$ensembl_gene_id),3]
   RnEG = as.character(HsENSG$RnEG)
   names(RnEG) = HsENSG[,1]
   RnEG = RnEG[!is.na(RnEG)]
   RnEG = RnEG[unique(names(RnEG))]
   # export data
   for (i in which(names(export) %in% names(RnEG))) {
    export[i] <- RnEG[names(RnEG) == names(export[i])] 
   }
   if(verbose) cat(length(RnEG))
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
   RnENSP = mget(HsENSP, hom.Hs.inpRATNO, ifnotfound = NA)
   RnENSP = sapply(RnENSP[!is.na(RnENSP)], function(x) {x[1]}) # select first match only
   if (length(RnENSP) == 0) stop("No RnENSP.")
   require("org.Rn.eg.db")
   RnEG = mget(RnENSP, org.Rn.egENSEMBLPROT2EG, ifnotfound = NA)
   RnEG = sapply(RnEG, function(x) {x[1]}) # select first match only
   RnEG = RnEG[!is.na(RnEG)]
   # match up Human to Rat results 
   names(RnEG) = names(HsENSP)[HsENSP %in% names(RnENSP)[RnENSP %in% names(RnEG)]]
   RnEG = RnEG[!is.na(names(RnEG))]
   # prepare data export
   export[names(export) %in% names(RnEG)] = RnEG
   if(verbose) cat(sum(names(export) %in% names(RnEG)))
   cat(".")
  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 }

 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 return(export)
}

