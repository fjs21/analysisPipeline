convertHumanToMouseGeneID <-
function(HsEG, verbose = F) {
 # Check for valid input
 if (length(HsEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 HsEG <- as.character(HsEG)

 cat("Converting human geneIds to mouse.")
 export = rep(0, length(HsEG))
 names(export) = HsEG  
 cat(".") 

 # Homologene based (Kuhn et al., 2008 - BMC)
 if (exists("homologene")) {
  require(annotationTools)
  suppressWarnings(MmEG <- getHOMOLOG(HsEG, 10090, homologene))
  MmEG <- sapply(MmEG, function(x) {as.character(x[1])}) 
  names(MmEG) = HsEG
  MmEG = MmEG[!is.na(MmEG)]
  export[names(export) %in% names(MmEG)] = MmEG
  if(verbose) cat(sum(names(export) %in% names(MmEG)))
  cat(".")
 } else {
  warning("Homologene is not loaded!\n")
 }

 # Ensembl BioMart based
 if (any(export == 0)) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Mm.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   HsENSG = getBM(attributes = c("entrezgene_id","ensembl_gene_id"), filters = "entrezgene_id",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   MmENSG = getBM(attributes = c("ensembl_gene_id","mmusculus_homolog_ensembl_gene"), filters = "entrezgene_id",
    values = HsEG[export == 0], mart = Hs.ensembl)
   if(dim(MmENSG)[1] == 0) stop("No MmENSG.")
   MmEG = getBM(attributes =  c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
    values = MmENSG[,2], mart = Mm.ensembl)
   MmEG = MmEG[!is.na(MmEG$entrezgene),]
   # match up Human to Mouse results
   MmENSG$MmEG = MmEG[match(MmENSG$mmusculus_homolog_ensembl_gene, MmEG$ensembl_gene_id),2]
   HsENSG$MmEG = MmENSG[match(HsENSG$ensembl_gene_id, MmENSG$ensembl_gene_id),3]
   MmEG = as.character(HsENSG$MmEG)
   names(MmEG) = HsENSG[,1]
   MmEG = MmEG[!is.na(MmEG)]
   MmEG = MmEG[unique(names(MmEG))]
   # export data
   for (i in which(names(export) %in% names(MmEG))) {
    export[i] <- MmEG[names(MmEG) == names(export[i])] 
   }
   if(verbose) cat(length(MmEG))
   cat(".")
  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 }

 # Bioconductor-based - depreciated
 # if (any(export == 0)) {
 #  require("org.Hs.eg.db")
 #  HsENSP = mget(HsEG[export == 0], org.Hs.egENSEMBLPROT, ifnotfound = NA)
 #  HsENSP = sapply(HsENSP[!is.na(HsENSP)], function(x) {x[1]}) # select first match only 
 #  require("hom.Hs.inp.db")
 #  tryCatch({
 #   if (length(HsENSP) == 0) stop("No HsENSP.")
 #   MmMGI = mget(HsENSP, hom.Hs.inpMUSMU, ifnotfound = NA)
 #   MmMGI = sapply(MmMGI[!is.na(MmMGI)], function(x) {x[1]}) # select first match only
 #   if (length(MmMGI) == 0) stop("No MmMGI.")   
 #   require("org.Mm.eg.db")
 #   MmEG = mget(MmMGI, org.Mm.egMGI2EG, ifnotfound = NA)
 #   MmEG = sapply(MmEG, function(x) {x[1]}) # select first match only
 #   MmEG = MmEG[!is.na(MmEG)]
 #   # match up Human to Mouse results 
 #   names(MmEG) = names(HsENSP)[HsENSP %in% names(MmMGI)[MmMGI %in% names(MmEG)]]
 #   MmEG = MmEG[!is.na(names(MmEG))]
 #   # prepare data export
 #   export[names(export) %in% names(MmEG)] = MmEG
 #   if(verbose) cat(sum(names(export) %in% names(MmEG)))
 #   cat(".")
 #  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 # }

 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 return(export)
}

