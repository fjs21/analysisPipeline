runAutomatedGSEA <-
function(selCoef, refSample, GSEAfolder = NULL, blnPrompt = FALSE, species = "Hs",
  pcutoff = 0.05, adjmethod = "fdr") {

 # function depends on contrast.matrix, eset
 if(!exists('eset') | !exists('contrast.matrix')) {
  stop("You must have a valid 'eset' and 'contrast.matrix' object before running this function",
   call. = FALSE)
 }
 
 # setup procedure
 if(is.null(GSEAfolder)) GSEAfolder = selCoef
    
 # set selCoef to numeric corresponding to colname of contrast.matrix (global variable)
 if (selCoef == "NULL") selCoef = NULL else {
  if(class(selCoef) != "numeric") {selCoef = which(colnames(contrast.matrix) == selCoef)} }

 # check refSample against eset (global variable)
 if(!any(eset$phenotype == refSample)) {
  setwd(wd)
  stop(paste("refSample '", refSample, "' does not exist in eset$phenotype.\n", sep=""), call. = FALSE)
 }
 if(length(which(eset$phenotype == refSample)) == 1) {
  setwd(wd)
  stop("GSEA does not support a single reference sample.", call. = FALSE)
 }

 if (!exists("chipAnnotation")) chipAnnotation = annotation(eset)
 require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
 require(GSEABase)
 require(genefilter)
 require(XML)

 # Filter eset to unique GeneIds
 if(!exists('nsF')) {
  nsF <<- nsFilter(eset, var.filter = FALSE)$eset # uses IQR as filter for duplicate probes
 }
 if(!exists('EG')) {
  # Retrieve probe names
  probes <<- featureNames(nsF)
  syms <<- as.character(unlist(mget(probes,eval(parse(text = paste(chipAnnotation,"SYMBOL",sep=""))))))
  EG <<- as.character(unlist(mget(probes,eval(parse(text = paste(chipAnnotation,"ENTREZID",sep=""))))))
 }

 if (species == "Mm") {
  EG <<- convertMouseToHumanGeneID(EG)
  names(probes) = names(EG)
  message(sum(EG != 0), " annotated with human GeneIDs, out of ", length(EG), " total.\n")
  # remove genes without homologs
  nsF <<- nsF[EG != 0,]
  probes <<- probes[EG != 0] # will still be original 
  #  syms <<- syms[EG != 0]
  EG <<- EG[EG != 0]
  
  # find possible duplicates
  testStat = genefilter:::rowIQRs(nsF)
  names(testStat) = names(EG)
  tSsp = split.default(testStat, EG)
  uniqEG <- sapply(tSsp, function(x) names(which.max(x)))
 
  nsF <<- nsF[names(EG) %in% uniqEG,]
  probes <<- probes[uniqEG]
  EG <<- EG[uniqEG]
     
  syms <<- unlist(mget(EG, org.Hs.egSYMBOL, ifnotfound = NA))
  # remove genes with no symbol
  EG <<- EG[which(!is.na(syms))]
  probes <<- probes[which(!is.na(syms))]
  nsF <<- nsF[which(!is.na(syms)),]
  syms <<- syms[which(!is.na(syms))]
 }

 # Not tested 06/08/12
 if (species == "Ss") {
  EG <<- convertPorcineToHumanGeneID(EG)
  message(sum(EG != 0), " annotated with human GeneIDs, out of ", length(EG), " total.\n")
  # remove genes without homologs
  nsF <<- nsF[EG != 0,]
  probes <<- probes[EG != 0] # will still be original 
#  syms <<- syms[EG != 0]
  EG <<- EG[EG != 0]
  syms <<- unlist(mget(EG, org.Hs.egSYMBOL, ifnotfound = NA))
 }

 # create and move to results directory
 wd = getwd()
 dir.create(paste("GSEA_",GSEAfolder,sep=""), showWarnings = FALSE)
 setwd(paste(wd,"/GSEA_",GSEAfolder,sep=""))
 
 # P-GSEA on Biocarta
 message("Biocarta GSEA analysis running...\n")
 data(bioCarta)

 featureNames(nsF) = EG
 GeneIDUniverse = unique(EG)

 gsl = list()
 for (i in 1:length(PathwayIDs)) {
  selected.Pathway = InteractionMap[InteractionIDs %in% PathwayMap[[i]]$Interactions]
  nodes.molecules = unique(as.character(unlist(sapply(selected.Pathway, function(n) c(n$inputs,n$agents,n$outputs,n$inhibitors)))))
  pathway.geneIDs = as.character(intersect(unique(unlist(sapply(MoleculeMap[nodes.molecules],function(n) n$GeneID))), GeneIDUniverse))
  gs = list(GeneSet(setName = PathwayNames[i], geneIds=pathway.geneIDs, geneIdType=EntrezIdentifier(),
    shortDescription = PathwayIDs[i]))
  gsl = c(gsl, gs)
 }
 gsc = GeneSetCollection(gsl)

 featureNames(nsF) = EG
 PGSEA.sig = run.GSEA(nsF = nsF, gsc = gsc, refSample = refSample, 
  file = "Biocarta_GSEA.PDF", blnPrompt = blnPrompt, range = c(5,500), coef = selCoef,
  pval = pcutoff, adjmethod = adjmethod)
 write.csv(PGSEA.sig, file = "Biocarta_GSEA.csv") 

 featureNames(nsF) = probes
 if(class(PGSEA.sig) != "character") {
  for (i in 1:min(length(PGSEA.sig$ID),25) ) {
   pathway = which(PathwayNames == PGSEA.sig$ID[i])
   print(PathwayNames[pathway])
   pdf(file=paste("Biocarta_", i,".pdf", sep=""), paper="USr", pointsize = 10, width=1600, height=900)
   bioCarta.plot(pathway,
    fit3$coefficients[rownames(fit3$coefficients) %in% probes,
    selCoef], species = species)
   bioCarta.heatmap(pathway, nsF, species = species)
   dev.off()
  }
 }

 # P-GSEA on Broad sets
 message("Broad C2-based GSEA analysis running...\n")
 # Retrieve GMT file

 if(exists('GSEA.C2')) {
  broad.gsc = getGmt(GSEA.C2, geneIdType=NullIdentifier(),
   collectionType=NullCollection(), sep="\t")

  featureNames(nsF)= syms
  PGSEA.sig = run.GSEA(nsF = nsF, gsc = broad.gsc, refSample = refSample, 
   file = "Broad_C2_GSEA.PDF", blnPrompt = blnPrompt,
   range = c(25,5000), coef = selCoef,
   pval = pcutoff, adjmethod = adjmethod) 

  if(class(PGSEA.sig) != "character") { 
   # get details of broad gene sets - define new class
   setClass("broadGSC", representation(NAME = "character", 
    PMID = "character", AUTHORS = "character",
    DESCRIPTION_BRIEF = "character", DESCRIPTION_FULL = "character"),
    prototype = prototype(PMID = "", AUTHORS = "",
    DESCRIPTION_BRIEF = "", DESCRIPTION_FULL = ""),
    where = topenv(parent.frame()) )

  library(RCurl)
   # on initialize use NAME to retrieve details
   setMethod("initialize", "broadGSC", function(.Object, ..., NAME) {
    tryCatch({
     broad.xml = getURL(asBroadUri(NAME,base="https://software.broadinstitute.org/gsea/msigdb/cards"))
     broad.xml = xmlTreeParse(broad.xml, useInternalNodes = T)
     broadGS = callNextMethod(.Object, ..., NAME = NAME, 
      PMID = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"PMID"),
      AUTHORS = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"AUTHORS"),
      DESCRIPTION_BRIEF = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"DESCRIPTION_BRIEF"),
      DESCRIPTION_FULL = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"DESCRIPTION_FULL"))
      #message("Parsed ",NAME,".",sep="")
      broadGS
      }, error = function(err) {
 	  message("Failed to parse XML for ", NAME, " pathway.")
	  .Object@NAME = NAME
	  .Object
    })
   }, where = topenv(parent.frame()) )
   # get details on significant pathways
   PGSEA.broadGSC = sapply(PGSEA.sig$ID, function(x) {
    new("broadGSC", NAME = x)})
   # add to table for export
   PGSEA.sig = data.frame(
    PMID = unlist(sapply(PGSEA.broadGSC, function(x) x@PMID)),
    AUTHORS = unlist(sapply(PGSEA.broadGSC, function(x) x@AUTHORS)),
    DESCRIPTION_BRIEF = unlist(sapply(PGSEA.broadGSC, function(x) x@DESCRIPTION_BRIEF)),
    PGSEA.sig)
  }
  write.csv(PGSEA.sig, file = "Broad_C2_GSEA.csv")
 } else {warning("Broad C2 file not found. Please specify location in 'GSEA.C2' global variable before running.")}
 
 # MOTIF examination C3
 message("Broad C3-based GSEA analysis running...\n")
 if (exists('GSEA.C3')) {
  broad.gsc = getGmt(GSEA.C3, geneIdType=NullIdentifier(),
   collectionType=NullCollection(), sep="\t")

  featureNames(nsF)= syms
  PGSEA.sig = run.GSEA(nsF = nsF, gsc = broad.gsc, refSample = refSample,
   file = "Broad_C3_GSEA.PDF", blnPrompt = blnPrompt,
   range = c(25,5000), coef = selCoef,
   pval = pcutoff, adjmethod = adjmethod)

  if(class(PGSEA.sig) != "character") {
   # get details of broad gene sets - define new class
   setClass("broadGSC", representation(NAME = "character", 
    PMID = "character", AUTHORS = "character",
    DESCRIPTION_BRIEF = "character", DESCRIPTION_FULL = "character"),
    prototype = prototype(PMID = "", AUTHORS = "",
    DESCRIPTION_BRIEF = "", DESCRIPTION_FULL = ""),
    where = topenv(parent.frame()) )
   # on initialize use NAME to retrieve details
   setMethod("initialize", "broadGSC", function(.Object, ..., NAME) {
    tryCatch({
      broad.xml = getURL(asBroadUri(NAME,base="https://software.broadinstitute.org/gsea/msigdb/cards"))
      broad.xml = xmlTreeParse(broad.xml, useInternalNodes = TRUE) 
     #message("Parsed ",NAME,".",sep="")
     callNextMethod(.Object, ..., NAME = NAME, 
      PMID = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"PMID"),
      AUTHORS = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"AUTHORS"),
      DESCRIPTION_BRIEF = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"DESCRIPTION_BRIEF"),
      DESCRIPTION_FULL = xmlGetAttr(getNodeSet(broad.xml,"//MSIGDB/GENESET")[[1]],"DESCRIPTION_FULL"))
    }, error = function(err) {
     message("Failed to parse XML for '", NAME, "' pathway.")
     .Object@NAME = NAME
     .Object
    })
   }, where = topenv(parent.frame()) )
  # get details on significant pathways
  PGSEA.broadGSC = sapply(PGSEA.sig$ID, function(x) {
   new("broadGSC", NAME = x)}) 
  # add to table for export
  PGSEA.sig = data.frame(
   PMID = unlist(sapply(PGSEA.broadGSC, function(x) x@PMID)),
   AUTHORS = unlist(sapply(PGSEA.broadGSC, function(x) x@AUTHORS)),
   DESCRIPTION_BRIEF = unlist(sapply(PGSEA.broadGSC, function(x) x@DESCRIPTION_BRIEF)),
    PGSEA.sig)
  }
  write.csv(PGSEA.sig, file = "Broad_C3_GSEA.csv")
 } else {warning("Broad C3 file not found. Please specify location in 'GSEA.C3' global variable before running.")}

 # KEGG and GO based GSEA
 message("KEGG GSEA analysis running...")
 data(KEGGGOgsc)
 
 featureNames(nsF)= EG
 PGSEA.sig = run.GSEA(nsF = nsF, gsc = KEGG.mcs, refSample = refSample,
  file = "KEGG_GSEA.PDF", blnPrompt = blnPrompt,
  range = c(5,1000), coef = selCoef,
  pval = pcutoff, adjmethod = adjmethod)
 write.csv(PGSEA.sig, file = "KEGG_GSEA.csv")

 if(class(PGSEA.sig) != "character") {
  featureNames(nsF)= probes
  for(i in 1:min(5,dim(PGSEA.sig)[1])) {
   KeggID = substr(PGSEA.sig$ID[i],1,regexpr("\\>",PGSEA.sig$ID[i])[1]-1)
   pdf(paste("KEGG_",KeggID,".pdf",sep=""), paper="USr", pointsize = 10, width=1600, height=900)
   tryCatch({
    pathway.probes = KEGG.plot(KeggID,
     fit3$coefficients[rownames(fit3$coefficients) %in% probes,selCoef], chipAnnotation, species = species)
    custom.heatmap(pathway.probes, nsF, title = PGSEA.sig$ID[i]) },
    error = function(e) cat("Error occured plotting kegg pathway '",KeggID,"'.") )
   dev.off()
  }
 }

 # GO GSEA
 message("GO GSEA analysis running...\n")
 featureNames(nsF)= EG
 PGSEA.sig = run.GSEA(nsF = nsF, gsc = GO.mcs, refSample = refSample,
  file = "GO_GSEA.PDF", range = c(5,1000), coef = selCoef,
  pval = pcutoff, adjmethod = adjmethod)
 write.csv(PGSEA.sig, file = "GO_GSEA.csv")

 # Clean-up
 #save.image("GSEA.RData")
 setwd(wd)
}

