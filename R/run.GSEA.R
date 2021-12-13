run.GSEA <-
function(nsF, gsc, refSample, file = "GSEA.pdf", 
 blnPrompt = FALSE,
 range = c(25,500), coef = NULL,
 skip = "NO", pval = 0.05, adjmethod = "BH"
) {

 require(PGSEA)

 maxLblWidth = max(nchar(as.character(nsF$phenotype)))
 #margin = c(1, 1, maxLblWidth/7*4, 20)	
 margin = c(1,1,3,10)

 # No longer used
 #pg <- PGSEA(exprs(nsF), gsc,
 # ref=which(nsF$phenotype == refSample), range=range)
 #plotRange = as.numeric(max(abs(range(pg, finite = TRUE))))
 #if(blnPrompt == TRUE) par(ask = TRUE)
 #smcPlot(pg[apply(pg, 1, function(x) any(!is.na(x))),], 
 # factor(nsF$phenotype), col = .rwb, scale = c(-plotRange,plotRange),
 # margin = margin, show.grid = TRUE, r.cex = 0.75, skip = skip)

 pgNF <- PGSEA(exprs(nsF), gsc, p.value = NA, 
  ref=which(nsF$phenotype == refSample), range = range)

 require(limma)
 design = get("design", envir=globalenv())
 if(exists('arrayw')) {
  PGSEA.fit <- lmFit(pgNF, design, weights = arrayw)
 } else {
  PGSEA.fit <- lmFit(pgNF, design)
 }
 PGSEA.fit <- contrasts.fit(PGSEA.fit, contrast.matrix)
 PGSEA.fit <- eBayes(PGSEA.fit)

# define 'desiredResults' matrix of desired results with rows 
# representing alternate contrast.matrix configurations
# and columns for each contrast in contrast.matrix
# 1 = UP regulated
# 0 = no change
# -1 = DOWN regulated
# NA = any regulation

# if 'desiredResults' not defined set to test for up and down regulation in selected
# coefficient only

 if(!exists("desiredResults")) { 
  desiredResults <- matrix(rep(NA,ncol(contrast.matrix) * 2), nrow = 2)
  desiredResults[1,coef] <- 1
  desiredResults[2,coef] <- -1
 } else {
  cat("Using 'desiredResults' matrix.\n")
 }

 dT = decideTests(PGSEA.fit, p.value = pval, adjust.method = adjmethod)
 testMatrix = apply(desiredResults, 1, function(dR) {
  apply(dT, 1, function(resultRow) { all(resultRow[!is.na(dR)] == dR[!is.na(dR)]) }) 
  })
 passfail = apply(testMatrix, 1, any)
 poi = rownames(as.matrix(dT[passfail,])) 
 PGSEA.sig <- topTable(PGSEA.fit, number=1000, p.value = 1, coef=coef)
 PGSEA.sig$ID <- rownames(PGSEA.sig)
 PGSEA.sig = PGSEA.sig[PGSEA.sig$ID %in% poi,]

 if(dim(PGSEA.sig)[1] == 0) {
  cat("No significant pathways.\n")
  return("no significant pathways")}
 if(dim(PGSEA.sig)[1] == 1) {
  cat("Only one significant pathway.\n")
  return(PGSEA.sig)}

 plotRange = as.numeric(max(abs(range(pgNF[PGSEA.sig$ID,], finite = TRUE))))

 smcPlot(pgNF[PGSEA.sig$ID,], 
  factor(nsF$phenotype), col = .rwb, scale = c(-plotRange,plotRange),
  margins = margin, show.grid = TRUE, r.cex = 0.75, skip = skip)
 par(ask = FALSE)

 pdf(file=file, paper="USr", pointsize = 10, width=1600, height=900)
 nsig = dim(PGSEA.sig)[1]
 cat(nsig, "significant pathways found.\n")
 plotRange = as.numeric(max(abs(range(pgNF[PGSEA.sig$ID,], finite = TRUE))))
 pageSize = 50
 pages = 1 + floor(nsig/pageSize)
 if(nsig%%pageSize == 1) {cat("Last pathway ommited from PDF output.")}
 for (page in 1:pages) {
  start = pageSize *(page-1) + 1
  end =  ifelse(page == pages, nsig%%pageSize + start - 1, pageSize * page)
  if(start < end) {
   smcPlot(pgNF[PGSEA.sig$ID[start:end],], 
    factor(nsF$phenotype), col = .rwb, scale = c(-plotRange,plotRange),
    margins = margin, show.grid = TRUE, r.cex = 0.75, skip = skip)
  }
 }
 dev.off()

 return(PGSEA.sig)
}

