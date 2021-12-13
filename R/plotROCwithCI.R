plotROCwithCI <- 
function (gs, fit3 = fit3, coef, main, pAUCthreshold = 0.2, add = FALSE, 
    ...) 
{
# plot confidence limits for single curve
    if (class(gs) == "GeneSet") {
        require(GSEABase)
        gs.EG <- geneIds(gs)
    }
    else {
        gs.EG <- as.character(gs)
    }
    if (!exists("EG") | !exists("probes")) {
        if (!exists("nsF")) {
            if (!exists("eset")) {
                stop("Please specify a valid eset object before continuing.")
            }
            require(genefilter)
            message("Generating geneID unique ExpressionSet...")
            nsF <<- nsFilter(eset, var.filter = FALSE)$eset
        }
        probes <<- featureNames(nsF)
        if (!exists("chipAnnotation")) 
            chipAnnotation <<- annotation(nsF)
        message("Annotating probes to unique geneID...")
        EG <<- getGeneIDFromProbeset(probes)
    }
    contrasts <- colnames(fit3$contrasts)
    data <- fit3$coefficients[probes, coef]
    truth <- (EG %in% gs.EG) * 1
    require(pROC)
    #rocobj <- roc(truth, data, direction = "<")
    rocobj <- plot.roc(truth, data, direction = "<", identity.col = 'red',
     lwd = 1, ci = TRUE, print.auc = TRUE,
     main = main, ...)
    abline(v = 1 - pAUCthreshold, col = "black", lty = 2)
   
    ciobj <- ci.se(rocobj,
     specificities=seq(0, 1, .05), boot.n = 100)
    plot(ciobj, type="shape", col="#1c61b6AA", lwd = 1)
    #plot(ci(rocobj, of="thresholds", thresholds="best"))     

    return(list(
     AUC = auc(rocobj), pAUC1 = auc(rocobj,
      partial.auc = c(1, 1 - pAUCthreshold)) ))
}