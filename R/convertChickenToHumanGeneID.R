convertChickenToHumanGeneID <-
function(GgEG, verbose = F) {
 # Check for valid input
 if (length(GgEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 GgEG <- as.character(GgEG)

 cat("Converting Chicken geneIds to human.")
 export = rep(0, length(GgEG))
 names(export) = GgEG  
 cat(".") 

 # Homologene based
 if (exists("homologene")) {
  require(annotationTools)
  suppressWarnings(HsEG <- getHOMOLOG(GgEG, 9606, homologene))
  HsEG <- sapply(HsEG, function(x) {as.character(x[1])}) 
  names(HsEG) = GgEG
  HsEG = HsEG[!is.na(HsEG)]
  export[names(export) %in% names(HsEG)] = HsEG
  if(verbose) cat(sum(names(export) %in% names(HsEG)))
  cat(".")
 } else {
  warning("Homologene is not loaded!\n")
 }

 # Ensembl BioMart based
 if (any(export == 0)) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Gg.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   GgENSG = getBM(attributes = c("entrezgene_id","ensembl_gene_id"), filters = "entrezgene_id",
    values = GgEG[export == 0], mart = Gg.ensembl)
   if(dim(GgENSG)[1] == 0) stop("No GgENSG.")
   HsENSG = getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), filters = "entrezgene_id",
    values = GgEG[export == 0], mart = Gg.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   HsEG = getBM(attributes =  c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
    values = HsENSG[,2], mart = Hs.ensembl)
   HsEG = HsEG[!is.na(HsEG$entrezgene),]
   # match up Human to Chicken results
   HsENSG$HsEG = HsEG[match(HsENSG$hsapiens_homolog_ensembl_gene, HsEG$ensembl_gene_id),2]
   GgENSG$HsEG = HsENSG[match(GgENSG$ensembl_gene_id, HsENSG$ensembl_gene_id),3]
   HsEG = as.character(GgENSG$HsEG)
   names(HsEG) = GgENSG[,1]
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


 # Bioconductor-based
 if (any(export == 0)) {
  require("org.Gg.eg.db")
  GgENSP = mget(GgEG[export == 0], org.Gg.egENSEMBLPROT, ifnotfound = NA)
  GgENSP = sapply(GgENSP[!is.na(GgENSP)], function(x) {x[1]}) # select first match only 
  require("hom.Hs.inp.db")
  tryCatch({
   if (length(GgENSP) == 0) stop("No GgENSP.")
   HsENSP = mget(GgENSP, revmap(hom.Hs.inpGALGA), ifnotfound = NA)
   HsENSP = sapply(HsENSP[!is.na(HsENSP)], function(x) {x[1]}) # select first match only
   if (length(HsENSP) == 0) stop("No HsENSP.")
   require("org.Hs.eg.db")
   HsEG = mget(HsENSP, org.Hs.egENSEMBLPROT2EG, ifnotfound = NA)
   HsEG = sapply(HsEG, function(x) {x[1]}) # select first match only
   HsEG = HsEG[!is.na(HsEG)]
   # match up Human to Chicken results 
   names(HsEG) = names(GgENSP)[GgENSP %in% names(HsENSP)[HsENSP %in% names(HsEG)]]
   HsEG = HsEG[!is.na(names(HsEG))]
   # prepare data export
   export[names(export) %in% names(HsEG)] = HsEG
   if(verbose) cat(sum(names(export) %in% names(HsEG)))
   cat(".")
  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 } 

 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 return(export)
}

