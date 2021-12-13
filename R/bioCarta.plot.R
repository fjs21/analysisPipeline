bioCarta.plot <-
function(pid, ratios = NULL, species = "Hs") { 

 if (!exists("InteractionMap")) data(bioCarta)

 #require(graph)
 if (!require(Rgraphviz)) {
  warning("Rgraphviz not available")
  return(0)
 }

 # Setup pathway - retrieve interactions & molecules
 selected.Pathway = InteractionMap[InteractionIDs %in% PathwayMap[[pid]]$Interactions]
 nodes.interactions = paste("i",names(selected.Pathway),sep="")
 nodes.molecules = unique(as.character(unlist(sapply(selected.Pathway, function(n) c(n$inputs,n$agents,n$outputs,n$inhibitors)))))

 # Interaction edges = outputs
 edges.interactions = list()
 for (i in 1:length(nodes.interactions)) {
	nodename.interactions = as.character( nodes.interactions[i] )
	edges.interactions[[nodename.interactions]] = list(edges = selected.Pathway[[i]]$outputs)
 }

 # Molecule edges = inputs, agents, inhibitors
 edges.molecules = list()
 for (i in 1:length(nodes.molecules)) {
	InteractionsWithMatchingInputNode =
	 which(sapply(selected.Pathway, function(n) any(n$inputs == nodes.molecules[i])|any(n$agents == nodes.molecules[i])|any(n$inhibitors == nodes.molecules[i]) ))

	if(length(InteractionsWithMatchingInputNode)>0) edges.molecules[[as.character(nodes.molecules[i])]] = list(edges = paste("i",names(selected.Pathway[InteractionsWithMatchingInputNode]), sep="")) 
	 else edges.molecules[[as.character(nodes.molecules[i])]] = list(edges = character(0))

 }
 # Generate graph object
 myNodes = c(nodes.interactions, nodes.molecules)
 myEdges = c(edges.interactions, edges.molecules)
 g <- new("graphNEL", nodes = myNodes, edgeL = myEdges, edgemode="directed") 

 # Node attributes
 # Generate labels & trim if necessary
 names(labels) = labels = nodes(g)
 labels = c(
   gsub(" ","\\\n",InteractionType[InteractionIDs %in% substr(nodes.interactions,2,nchar(nodes.interactions))])
 ,
   as.character(sapply(MoleculeMap[nodes.molecules], function(n) {
	ifelse(length(n$Symbol)==0,paste(substr( gsub("/","\\\n",n$Alias),1,30), collapse = "\n"), n$Symbol)
	})))
 names(labels) = nodes(g)
 nc = nchar(labels)
 # other node attributes
 fill = c(rep("white", times=length(nodes.interactions)),rep("lightblue", times=length(nodes.molecules)))
 shape = c(rep("box", times=length(nodes.interactions)),rep("ellipse", times=length(nodes.molecules)))
 fixedsize = c(rep(TRUE, times=length(nodes.interactions)),rep(FALSE, times=length(nodes.molecules)))
 width = c(rep(20, times=length(nodes.interactions)),rep(75, times=length(nodes.molecules)))
 height = c(rep(20, times=length(nodes.interactions)),rep(50, times=length(nodes.molecules)))
 fontscale = 12
 fontsize = fontscale * max(nchar(unlist(strsplit(labels,"\n", fixed=TRUE)))) /
  sapply(labels, function(x) {
   max(nchar(unlist(strsplit(x,"\n", fixed=TRUE))))
 })
 fontsize[is.infinite(fontsize)] = fontscale
 names(fill) = names(shape) = names(fontsize) = names(width) = names(height) = nodes(g)

 # Edge Antributes
 arrowsize = c(rep(1, times=length(edgeNames(g))))
 names(arrowsize) = edgeNames(g)
 edgeColor = NULL
 # Molecule edges = agent & output
 for (i in 1:length(nodes.molecules)) {
	InteractionsWithMatchingInputNode = which(sapply(selected.Pathway, function(n) {
	 any(n$agents == nodes.molecules[i]) | any(n$outputs == nodes.molecules[i])}
       ))
	if(length(InteractionsWithMatchingInputNode)>0) {
		agentedges = paste(nodes.molecules[i],names(InteractionsWithMatchingInputNode),sep="~i")
		color.temp = rep("green", times=length(agentedges))
		names(color.temp) = agentedges
		edgeColor = c(edgeColor,color.temp)
		arrowhead.temp = rep("normal", times=length(agentedges))
		names(arrowhead.temp) = agentedges
		arrowhead = c(arrowhead,arrowhead.temp)
	}
 }
 # Molecule edges = inhibitor
 for (i in 1:length(nodes.molecules)) {
	InteractionsWithMatchingInputNode = which(sapply(selected.Pathway, function(n) any(n$inhibitors == nodes.molecules[i])))
	if(length(InteractionsWithMatchingInputNode)>0) {
		inhibitoredges = paste(nodes.molecules[i],names(InteractionsWithMatchingInputNode),sep="~i")
		color.temp = rep("red", times=length(inhibitoredges))
		names(color.temp) = inhibitoredges
		edgeColor = c(edgeColor,color.temp)
		arrowhead.temp = rep("tee", times=length(inhibitoredges))
		names(arrowhead.temp) = inhibitoredges
		arrowhead = c(arrowhead,arrowhead.temp)
	}
 }

if(!is.null(ratios)) { 
 # Calculate node colors according to "ratios"
 require(AnnotationDbi)
 require(paste(chipAnnotation,".db",sep=""), character.only = TRUE)
 pathway.geneIDs = sapply(MoleculeMap[nodes.molecules],function(n) n$GeneID)

 if (species == "Mm") {
  pathway.geneIDs <- sapply(pathway.geneIDs, function(x) {
   x = convertHumanToMouseGeneID(x)
   x = x[x != 0] 
  })
 }

 pathway.geneIDs = pathway.geneIDs[!is.na(pathway.geneIDs)]

 # EQUAL WEIGHTING PROBE SETS & GENES - disabled
 # pathway.probes = sapply(pathway.geneIDs, function(x) { 
 #  as.character(unlist(mget(x, revmap(eval(parse(text = paste(chipAnnotation,"ENTREZID",sep="")))), ifnotfound = NA)))})
 # pathway.probes = pathway.probes[!is.na(pathway.probes)]

 # Only show PIE with multiple genes per node - otherwise mean
 pathway.probes = sapply(pathway.geneIDs, function(x) { 
  mget(x, revmap(eval(parse(text = paste(chipAnnotation,"ENTREZID",sep="")))), ifnotfound = NA)})

 # pathway.expression_new = sapply(pathway.probes, function(x) {
 #   ratios[intersect(as.character(x), names(ratios))] })

 pathway.expression_new = sapply(pathway.probes,function(x) {
  if (length(x) == 0) return(0)
  y = sapply(x, function(x) {
   mean(ratios[intersect(as.character(x), names(ratios))])})
  y = y[!is.nan(y)]
  return(y)
 })

 pathway.expression = sapply(pathway.expression_new, function(x) {
   mean(x)})
 pathway.expression = pathway.expression[!is.na(pathway.expression)]

 hmcol <- colorRampPalette(brewer.pal(10, "RdYlGn")[10:1])(256)  

 pathway.fill = c("darkgreen","green","palegreen","lightyellow","orange","red","darkred")[cut(pathway.expression,
  breaks = c(-Inf, -2, -1, -0.5, 0.5, 1, 2, +Inf))]
 names(pathway.fill) = names(pathway.expression)
 fill[names(pathway.fill)] = pathway.fill
}

g <- layoutGraph(g, layoutType = "dot")
nodeRenderInfo(g) = list(label = labels, shape = shape, fontsize = fontsize,
 height = height, lWidth = width/2, rWidth = width/2,
 fill = fill)
# enable pie glyphs
expressionNode <- function(x, label, lWidth, rWidth, height, ...) {
 expressionData = pathway.expression_new[[names(label)]]
 fill = c("darkgreen","green","palegreen","lightyellow","orange","red","darkred")[cut(expressionData,
  breaks = c(-Inf, -2, -1, -0.5, 0.5, 1, 2, +Inf))]
 nodeDiameter = min( diff(x[,1]), diff(x[,2]) ) 
 pieGlyph (rep(1,times = length(expressionData)),
  xpos = mean(x[,1]), ypos = mean(x[,2]),
  radius = nodeDiameter/2, col = fill, labels = NA )
}

if(!is.null(ratios)) { 
 nodesWithPie = names(pathway.expression_new)[sapply(pathway.expression_new, length) > 1]
 nodeRenderInfo(g)$shape[nodesWithPie] = sapply(nodeRenderInfo(g)$shape[nodesWithPie], function(x) expressionNode)
}

edgeRenderInfo(g) = list(col = edgeColor, arrowhead = arrowhead) # arrowsize = arrowsize, 
graphRenderInfo(g) = list(graph=list(rankdir="TB"), main = PathwayNames[pid])
g <- renderGraph(g)
}

