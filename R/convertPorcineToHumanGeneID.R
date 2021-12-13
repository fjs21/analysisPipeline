convertPorcineToHumanGeneID <-
function(SsEG, verbose = F) {
 # Check for valid input
 if (length(SsEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 SsEG <- as.character(SsEG)

 cat("Converting porcine geneIds to human.")
 export = rep(0, length(SsEG))
 names(export) = SsEG  
 cat(".") 

# Homologene based - no Sus Scrofa records (Tax ID:9823 in homologene table)
 if (exists("homologene")) {
  require(annotationTools) 
  suppressWarnings(HsEG <- getHOMOLOG(SsEG, 9606, homologene))
  HsEG <- sapply(HsEG, function(x) {as.character(x[1])}) 
  names(HsEG) = SsEG
  HsEG = HsEG[!is.na(HsEG)]
  export[names(export) %in% names(HsEG)] = HsEG
  if(verbose) cat(sum(names(export) %in% names(HsEG)))
  cat(".")
 } else {
  warning("Homologene is not loaded!")
 }

# Ensembl BioMart based
 if (any(any(export == 0))) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Ss.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   SsENSG = getBM(attributes = c("entrezgene_id","ensembl_gene_id"), filters = "entrezgene_id",
    values = SsEG[export == 0], mart = Ss.ensembl)
   if(dim(SsENSG)[1] == 0) stop("No SsENSG.")
   HsENSG = getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), filters = "entrezgene_id",
    values = SsEG[export == 0], mart = Ss.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   HsEG = getBM(attributes =  c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
    values = HsENSG[,2], mart = Hs.ensembl)
   HsEG = HsEG[!is.na(HsEG$entrezgene),]
   # match up Human to Mouse results
   HsENSG$HsEG = HsEG[match(HsENSG$hsapiens_homolog_ensembl_gene, HsEG$ensembl_gene_id),2]
   SsENSG$HsEG = HsENSG[match(SsENSG$ensembl_gene_id, HsENSG$ensembl_gene_id),3]
   HsEG = as.character(SsENSG$HsEG)
   names(HsEG) = SsENSG[,1]
   HsEG = HsEG[!is.na(HsEG)]
   HsEG = HsEG[unique(names(HsEG))]
   # export data
   for (i in which(names(export) %in% names(HsEG))) {
    export[i] <- HsEG[names(HsEG) == names(export[i])] 
   }
   if(verbose) cat(length(HsEG))
   cat(".")
  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 }

# Bioconductor-based - does not work 6/8/12 as no SUSSC records available in hom.Hs.inp
# if (any(any(export == 0))) {
#  require("org.Ss.eg.db")
#  SsENSP = mget(SsEG[export == 0], org.Ss.egENSEMBLPROT, ifnotfound = NA)
#  SsENSP = sapply(SsENSP[!is.na(SsENSP)], function(x) {x[1]}) # select first match only 
#  require("hom.Hs.inp.db")
#  tryCatch({
#   if (length(MmENSP) == 0) stop("No SsENSP.")
#   HsENSP = mget(SsENSP, revmap(hom.Hs.inpSUSSC), ifnotfound = NA)
#   HsENSP = sapply(HsENSP[!is.na(HsENSP)], function(x) {x[1]}) # select first match only
#   if (length(HsENSP) == 0) stop("No HsENSP.")
#   require("org.Hs.eg.db")
#   HsEG = mget(HsENSP, org.Hs.egENSEMBLPROT2EG, ifnotfound = NA)
#   HsEG = sapply(HsEG, function(x) {x[1]}) # select first match only
#   HsEG = HsEG[!is.na(HsEG)]
#   # match up Human to Porcine results 
#   names(HsEG) = names(SsENSP)[SsENSP %in% names(HsENSP)[HsENSP %in% names(HsEG)]]
#   HsEG = HsEG[!is.na(names(HsEG))]
#   # prepare data export
#   export[names(export) %in% names(HsEG)] = HsEG
#   if(verbose) cat(sum(names(export) %in% names(HsEG)))
#   cat(".")
#  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
# }

 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 return(export)
}

