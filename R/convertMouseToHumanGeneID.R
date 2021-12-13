convertMouseToHumanGeneID <-
function(MmEG, verbose = F) {
 # Check for valid input
 if (length(MmEG) == 0) {
  return(character(0))
 }
 if(!exists('Hs.ensembl')) setupFindHomologs()
 MmEG <- as.character(MmEG)

 cat("Converting mouse geneIds to human.")
 export = rep(0, length(MmEG))
 names(export) = MmEG  
 cat(".") 

 # Homologene based
 if (exists("homologene")) {
   tryCatch({
    require(annotationTools)
    suppressWarnings(HsEG <- getHOMOLOG(MmEG, 9606, homologene))
    HsEG <- unlist(sapply(HsEG, function(x) {
      if(length(x)==1) return(as.character(x[1])) else return(NA)
      }) )
    #if (is.null(HsEG)) stop("Homologene: HsEG is null.")
    names(HsEG) = MmEG
    HsEG = HsEG[!is.na(HsEG)]
    HsEG = HsEG[isUnique(HsEG)]
    export[names(export) %in% names(HsEG)] = HsEG
    if(verbose) cat(" Hom: ",sum(names(export) %in% names(HsEG)))
    cat(".")
   }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } ) 
 } else {
  warning("Homologene is not loaded!\n")
 }

 # Ensembl BioMart based
 if (any(export == 0)) {
  require(biomaRt)
  if (!exists("Hs.ensembl")|!exists("Mm.ensembl")) stop("Biomart is not loaded")
  tryCatch({
   MmENSG = getBM(attributes = c("entrezgene_id","ensembl_gene_id"), filters = "entrezgene_id",
    values = MmEG[export == 0], mart = Mm.ensembl)
   MmENSG = MmENSG[isUnique(MmENSG$ensembl_gene_id),]
   if(dim(MmENSG)[1] == 0) stop("No MmENSG.")
   HsENSG = getBM(attributes = c("ensembl_gene_id","hsapiens_homolog_ensembl_gene"), filters = "entrezgene_id",
    values = MmEG[export == 0], mart = Mm.ensembl)
   if(dim(HsENSG)[1] == 0) stop("No HsENSG.")
   HsEG = getBM(attributes =  c("ensembl_gene_id","entrezgene_id"), filters = "ensembl_gene_id",
    values = HsENSG[,2], mart = Hs.ensembl)
   HsEG = HsEG[!is.na(HsEG$entrezgene),]
   # match up Human to Mouse results
   HsENSG$HsEG = HsEG[match(HsENSG$hsapiens_homolog_ensembl_gene, HsEG$ensembl_gene_id),2]
   # remove one-to-multi orthologs in Ensembl
   HsENSG = HsENSG[!HsENSG$HsEG %in% HsENSG$HsEG[duplicated(HsENSG$HsEG)],]
   if (dim(HsENSG)[1] == 0) 
     stop("No HsENSG after removing one-to-multi orthologs.")
   MmENSG$HsEG = HsENSG[match(MmENSG$ensembl_gene_id, HsENSG$ensembl_gene_id),3]
   HsEG = as.character(MmENSG$HsEG)
   names(HsEG) = MmENSG[,1]
   HsEG = HsEG[!is.na(HsEG)]
   HsEG = HsEG[unique(names(HsEG))]
   # remove would be duplicate annotations
   HsEG = HsEG[isUnique(HsEG)]
   HsEG = HsEG[!HsEG %in% export]
   # export data
   for (i in which(names(export) %in% names(HsEG))) {
    export[i] <- HsEG[names(HsEG) == names(export[i])] 
   }
   if(verbose) cat(" Ens: ",length(HsEG))
   cat(".")
  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 }

 # Bioconductor-based - depreciated
 # if (any(export == 0)) {
 #  require("org.Mm.eg.db")
 #  MmENSP = mget(MmEG[export == 0], org.Mm.egENSEMBLPROT, ifnotfound = NA)
 #  MmENSP = sapply(MmENSP[!is.na(MmENSP)], function(x) {x[1]}) # select first match only 
 #  require("hom.Mm.inp.db")
 #  tryCatch({
 #   if (length(MmENSP) == 0) stop("No MmENSP.")
 #   HsENSP = mget(MmENSP, hom.Mm.inpHOMSA, ifnotfound = NA)
 #   HsENSP = sapply(HsENSP[!is.na(HsENSP)], function(x) {x[1]}) # select first match only
 #   HsENSP 
 #   if (length(HsENSP) == 0) stop("No new BioC homologs found.")
 #   require("org.Hs.eg.db")
 #   HsEG = mget(HsENSP, org.Hs.egENSEMBLPROT2EG, ifnotfound = NA)
 #   HsEG = sapply(HsEG, function(x) {x[1]}) # select first match only
 #   HsEG = HsEG[!is.na(HsEG)]
 #   HsEG = HsEG[isUnique(HsEG)]
 #   # match up Human to Mouse results 
 #   names(HsEG) = names(MmENSP)[MmENSP %in% names(HsENSP)[HsENSP %in% names(HsEG)]]
 #   HsEG = HsEG[!is.na(names(HsEG))]
 #   # remove would be duplicate annotations
 #   HsEG = HsEG[!HsEG %in% export]
 #   # prepare data export
 #   export[names(export) %in% names(HsEG)] = HsEG
 #   if(verbose) cat(" bioC: ",sum(names(export) %in% names(HsEG)))
 #   cat(".")
 #  }, error = function(e) { if(verbose) cat(conditionMessage(e)) else cat("e") } )
 # } 

 cat("done.", ceiling(sum(export != 0)/length(export)*100) , "% annotated.\n")
 if(verbose) if(anyDuplicated(export[export != 0])) warning("Duplicate annotations!")
 return(export)
}

