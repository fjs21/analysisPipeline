bioCarta.heatmap <-
function (pwNumber, nsF, species = "Hs") {

 if (!is(nsF, "ExpressionSet")) {
  stop("bioCarta.heatmap only supports Expression set objects in 'nsF'.")
 }
 if (!exists("PathwayNames")) data(bioCarta)
 if (!exists("chipAnnotation")) assign("chipAnnotation", annotation(eset), .GlobalEnv)

 require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
 require("bioDist")

 pwName = PathwayNames[pwNumber]
 selected.Pathway = InteractionMap[InteractionIDs %in% PathwayMap[[pwNumber]]$Interactions]
 nodes.interactions = paste("i",names(selected.Pathway),sep="")
 nodes.molecules = unique(as.character(unlist(sapply(selected.Pathway, function(n) c(n$inputs,n$agents,n$outputs,n$inhibitors)))))

 pathway.geneIDs = unlist(sapply(MoleculeMap[nodes.molecules],function(n) n$GeneID))
 
 if (species == "Mm") {
  pathway.geneIDs = convertHumanToMouseGeneID(pathway.geneIDs)
  pathway.geneIDs = pathway.geneIDs[pathway.geneIDs != 0]
 }

 pathway.probes = as.character(unlist(sapply(pathway.geneIDs, function(x)
  tryCatch(get(x, revmap(eval(parse(text = paste(chipAnnotation,"ENTREZID",sep=""))))), error = function(n) print(n)) )))

 custom.heatmap(pathway.probes, nsF, paste("Heatmap - \n", pwName))
}

