plotROC <-
 function(gs, fit3 = fit3, coef, main, pAUCthreshold = 0.2, add = FALSE, 
  depletion = FALSE, ...) {

 # if 'gs' is a GeneSet recover geneIds
 if (class(gs) == "GeneSet") {
  require(GSEABase)
  gs.EG <- geneIds(gs)
 } else {
  gs.EG <- as.character(gs)
 } 

 # determine if global variables 'EG' and 'probes' objects exist 
 # if not, attempt to generate from 1st 'nsF' and then 'eset' objects
 if (!exists("EG") | !exists("probes")) {
  if (!exists("nsF")) {
   if (!exists("eset")) {stop("Please specify a valid eset object before continuing.")}
   require(genefilter)
   message("Generating geneID unique ExpressionSet...")
   nsF <<- nsFilter(eset, var.filter = FALSE)$eset
  }
  probes <<- featureNames(nsF)
  if(!exists("chipAnnotation")) chipAnnotation <<- annotation(nsF)
  message("Annotating probes to unique geneID...")
  EG <<- getGeneIDFromProbeset(probes)  
 }

 truth <- (EG %in% gs.EG) * 1
 data <- fit3$coefficients[probes,coef]
 if (depletion) data = data * -1
 require(ROC)
 roc <- rocdemo.sca(truth, data,
  markerLabel = "mL", caseLabel = "cL")
 if(add == FALSE) {
  plot(roc, line=TRUE, ylab = "sensitivity", xlab = "(1 - specificity)", ...)
  lines(c(0,1),c(0,1), col = "red")
  abline(v=pAUCthreshold, col = "black", lty = 2) 
  legend("bottomright",paste("  AUC = ",round(ROC::AUC(roc),digits = 2),
   "\npAUC = ", round(ROC::pAUC(roc, pAUCthreshold), digits = 3)), box.lty = 0, bty = 'n')
  title(main)
 } else {
  plot(roc, line=TRUE, add = T, ...)
 }

 return(list(AUC = ROC::AUC(roc), pAUC = ROC::pAUC(roc, pAUCthreshold),
  truth = truth, data = data))
}

