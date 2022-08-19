KEGG.plot <-
  function(KeggID,
           ratios = NULL,
           chipAnnotation = "hgu133plus2",
           expandGenes = FALSE,
           species = "Hs") {
    require(KEGGgraph)
    if (!require(Rgraphviz)) {
      warning("Rgraphviz not available")
      return(0)
    }
    require(org.Hs.eg.db)
    require(paste(chipAnnotation, ".db", sep = ""), character.only = TRUE)
    
    #KEGGurl = gsub("_","-",getKGMLurl(KeggID, organism = "hsa"))
    #KEGGurl = getCategoryIndepKGMLurl(KeggID, organism = "hsa")
    # Update July 19, 2011 - new KEGGurl avoids FTP site
    KEGGurl = paste(
      "http://www.genome.jp/kegg-bin/download?entry=hsa",
      KeggID,
      "&format=kgml",
      sep = ""
    )
    KEGGurl = paste("https://rest.kegg.jp/get/hsa",
                    KeggID, "/kgml", sep = "")
    download.file(KEGGurl, "KGML.xml")
    KEGGpathway = parseKGML("KGML.xml")
    keggG = KEGGpathway2Graph(KEGGpathway, genesOnly = FALSE, expandGenes = expandGenes)
    
    # remove isolated nodes
    outs <- sapply(edges(keggG), length) > 0
    ins <- sapply(inEdges(keggG), length) > 0
    ios <- outs | ins
    keggG <- subGraph(nodes(keggG)[ios], keggG)
    if (length(nodes(keggG)) == 0) {
      keggG = KEGGpathway2Graph(KEGGpathway,
                                genesOnly = FALSE,
                                expandGenes = TRUE)
      pathway.geneIDs <- translateKEGG2GeneID(nodes(keggG))
      if (species == "Mm") {
        pathway.geneIDs <- convertHumanToMouseGeneID(pathway.geneIDs)
      }
      pathway.probes = sapply(pathway.geneIDs, function(x) {
        mget(x, revmap(eval(parse(
          text = paste(chipAnnotation, "ENTREZID", sep = "")
        ))), ifnotfound = NA)
      })
      names(pathway.probes) = nodes(keggG)
      return(unlist(pathway.probes))
    }
    
    keggGnodedata <- getKEGGnodeData(keggG)
    keggGedgedata <- getKEGGedgeData(keggG)
    
    # get Edge labels (current only phosphorylation supported)
    edgeTypes = sapply(keggGedgedata, function(x) {
      paste(unlist(sapply(getSubtype(x), function(y)
        getName(y))), collapse = "_")
    })
    edgeLabel = ifelse(regexpr("phosphorylation", edgeTypes) > 0,
                       "P*",
                       as.character(""))
    names(edgeLabel) = names(edgeTypes)
    
    # layout graph
    keggG <- layoutGraph(keggG,
                         edgeAttrs = list(label = edgeLabel),
                         layoutType = "dot")
    
    ## translate the KEGG IDs into Gene Symbol
    if (expandGenes) {
      pathway.geneIDs <- translateKEGGID2GeneID(nodes(keggG))
      if (species == "Mm") {
        pathway.geneIDs <- convertHumanToMouseGeneID(pathway.geneIDs)
      }
      nodesNames <-
        sapply(mget(pathway.geneIDs, org.Hs.egSYMBOL, ifnotfound = NA),
               "[[",
               1)
      nodesNames[is.na(as.numeric(pathway.geneIDs))] =
        sapply(keggGnodedata[nodes(keggG)[is.na(as.numeric(pathway.geneIDs))]],
               function(x)
                 getDisplayName(x))
    } else {
      pathway.geneIDs <- sapply(keggGnodedata, function(x) {
        translateKEGGID2GeneID(x@name)
      })[nodes(keggG)]
      if (species == "Mm") {
        pathway.geneIDs <- sapply(pathway.geneIDs, function(x) {
          x = convertHumanToMouseGeneID(x)
          x = x[x != 0]
        })
      }
      #  nodesNames <- unlist(sapply(keggGnodedata, function(x) {
      #   getDisplayName(x) }))[nodes(keggG)]
      nodesNames <- getDisplayName(keggG)
    }
    
    # set up Node Attributes
    nodesNames = gsub(" ", "\\\n", nodesNames)
    if (expandGenes) {
      nodesFill = ifelse(is.na(as.numeric(pathway.geneIDs)), "purple", "lightblue")
      nodesShape = ifelse(is.na(as.numeric(pathway.geneIDs)), "rectangle", "ellipse")
    } else {
      nodesFill = sapply(pathway.geneIDs, function(x) {
        ifelse(all(is.na(as.numeric(x))), "purple", "lightblue")
      })
      nodesShape = sapply(pathway.geneIDs, function(x) {
        ifelse(all(is.na(as.numeric(x))), "rectangle", "ellipse")
      })
    }
    # Take default font size (14) and scale by fraction of length vs. median
    # multiply by number of line breaks introduced
    nodesFontSize = 24 * median(nchar(nodesNames)) / nchar(nodesNames) *
      sapply(gregexpr("\\\n", nodesNames), function(x)
        ifelse(any(x == -1), 1, length(x) + 1))
    names(nodesFontSize) = names(nodesShape) = names(nodesFill) = names(nodesNames) = nodes(keggG)
    
    # Change color if ratio information available
    if (!is.null(ratios)) {
      # EQUAL WEIGHTING PROBE SETS & GENES - disabled
      #  pathway.probes = sapply(pathway.geneIDs, function(x) {
      #   as.character(unlist(mget(x, revmap(eval(parse(text = paste(chipAnnotation,"ENTREZID",sep="")))), ifnotfound = NA)))})
      #  names(pathway.probes) = nodes(keggG)
      #  pathway.probes = pathway.probes[!is.na(pathway.probes)]
      #  pathway.expression_new = sapply(pathway.probes, function(x) {
      #   ratios[intersect(as.character(x), names(ratios))] })
      
      # ONE PIE SEGMENT PER GENE (PROBE SETS AVERAGED)
      pathway.probes = sapply(pathway.geneIDs, function(x) {
        mget(x, revmap(eval(parse(
          text = paste(chipAnnotation, "ENTREZID", sep = "")
        ))), ifnotfound = NA)
      })
      names(pathway.probes) = nodes(keggG)
      pathway.expression_new = sapply(pathway.probes, function(x) {
        y = sapply(x, function(x) {
          mean(ratios[intersect(as.character(x), names(ratios))])
        })
        y = y[!is.nan(y)]
      })
      pathway.expression = sapply(pathway.expression_new, function(x) {
        mean(x)
      })
      pathway.expression = pathway.expression[!is.na(pathway.expression)]
      
      pathway.fill = c("darkgreen",
                       "green",
                       "palegreen",
                       "lightyellow",
                       "orange",
                       "red",
                       "darkred")[cut(pathway.expression,
                                      breaks = c(-Inf,-2,-1,-0.5, 0.5, 1, 2,+Inf))]
      names(pathway.fill) = names(pathway.expression)
      nodesFill[names(pathway.fill)] = pathway.fill
    }
    
    nodeRenderInfo(keggG) <-
      list(
        label = nodesNames,
        shape = nodesShape,
        fill = nodesFill,
        fontsize = nodesFontSize
      )
    
    # enable pie glyphs
    expressionNode <- function(x, label, lWidth, rWidth, height, ...) {
      expressionData = pathway.expression_new[[names(label)]]
      fill = c("darkgreen",
               "green",
               "palegreen",
               "lightyellow",
               "orange",
               "red",
               "darkred")[cut(expressionData,
                              breaks = c(-Inf,-2,-1,-0.5, 0.5, 1, 2,+Inf))]
      nodeDiameter = min(diff(x[, 1]), diff(x[, 2]))
      pieGlyph (
        rep(1, times = length(expressionData)),
        xpos = mean(x[, 1]),
        ypos = mean(x[, 2]),
        radius = nodeDiameter / 2,
        col = fill,
        labels = NA
      )
    }
    
    if (!is.null(ratios)) {
      nodesWithPie = names(pathway.expression_new)[sapply(pathway.expression_new, length) > 1]
      nodeRenderInfo(keggG)$shape[nodesWithPie] = sapply(nodeRenderInfo(keggG)$shape[nodesWithPie],
                                                         function(x)
                                                           expressionNode)
    }
    
    edgeColor = factor(edgeTypes)
    edgeArrowhead = factor(edgeTypes)
    colorKey = c(
      "green",
      "palegreen",
      "green",
      "yellow",
      "purple",
      "grey",
      "cyan",
      "palegreen",
      "grey",
      "red",
      "red",
      "grey",
      "grey",
      "orange",
      "orange",
      "red",
      "yellow",
      "orange"
    )
    arrowheadKey = c(
      "normal",
      "normal",
      "normal",
      "dot",
      "diamond",
      "normal",
      "odot",
      "normal",
      "normal",
      "tee",
      "tee",
      "normal",
      "none",
      "diamond",
      "diamond",
      "tee",
      "dot",
      "diamond"
    )
    #### to do: add dashed lines for indirect ####
    names(colorKey) = names(arrowheadKey) = c(
      "activation",
      "activation_indirect",
      "activation_phosphorylation",
      "binding/association",
      "compound",
      "dephosphorylation",
      "dissociation",
      "expression",
      "indirect",
      "inhibition",
      "inhibition_phosphorylation",
      "phosphorylation",
      "missing",
      "ubiquitination",
      "activation_ubiquination",
      "inhibition_dephosphorylation",
      "phosphorylation_binding/association",
      "ubiquination_inhibition"
    )
    #levels(edgeColor) = colorKey
    #levels(edgeArrowhead) = arrowheadKey
    edgeColor = colorKey[match(edgeColor, names(colorKey))]
    edgeArrowhead = arrowheadKey[match(edgeArrowhead, names(arrowheadKey))]
    names(edgeColor) = names(edgeArrowhead) = names(edgeTypes)
    
    edgeRenderInfo(keggG) <- list(
      arrowsize = "0.7",
      col = edgeColor,
      fontsize = 10,
      textCol = edgeColor,
      arrowhead = edgeArrowhead
    )
    
    graphRenderInfo(keggG) = list(graph = list(rankdir = "TB"),
                                  main = getTitle(KEGGpathway))
    renderGraph(keggG)
    
    if (exists('pathway.probes'))
      return(unlist(pathway.probes))
  }
