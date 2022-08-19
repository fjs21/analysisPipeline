setupFindHomologs <-
function() {
 require(annotationTools)
 require(biomaRt)
 message("Setting up annotation functions.")
 if (exists('homologeneFile')) { 
  homologene <<- read.delim(homologeneFile, header = FALSE)
 } else {
  warning("Please specify a homologeneFile if you wish to take advantage of NCBI homologene annotations.")
 }
 Hs.ensembl <<- useMart("ENSEMBL_MART_ENSEMBL",dataset="hsapiens_gene_ensembl",host="https://www.ensembl.org")
 #useMart("ensembl", dataset="hsapiens_gene_ensembl")
 Mm.ensembl <<- useMart("ENSEMBL_MART_ENSEMBL",dataset="mmusculus_gene_ensembl",host="https://www.ensembl.org")
 Rn.ensembl <<- useMart("ENSEMBL_MART_ENSEMBL",dataset="rnorvegicus_gene_ensembl",host="https://www.ensembl.org")
 Bt.ensembl <<- useMart("ENSEMBL_MART_ENSEMBL",dataset="btaurus_gene_ensembl",host="https://www.ensembl.org")
 Ss.ensembl <<- useMart("ENSEMBL_MART_ENSEMBL",dataset="sscrofa_gene_ensembl",host="https://www.ensembl.org")
 Gg.ensembl <<- useMart("ENSEMBL_MART_ENSEMBL",dataset="ggallus_gene_ensembl",host="https://www.ensembl.org")
 message("Done.")
}

