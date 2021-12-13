ROC.compare2coef <- 
function (gs, fit3 = fit3, coef1, coef2, main, pAUCthreshold = 0.2, add = FALSE, 
    ...) 
{
    require(verification)
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
    truth <- (EG %in% gs.EG) * 1
    data1 <- fit3$coefficients[probes, coef1]
    require(pROC)
    rocobj1 <- roc(truth, data1, direction = "<")
    plot.roc(rocobj1, col = "#1c61b6", identity.col = 'red',
     lwd = 1, main = main, ...)
    abline(v = 1 - pAUCthreshold, col = "black", lty = 2)
    data2 <- fit3$coefficients[probes, coef2]
    rocobj2 <- roc(truth, data2, direction = "<")
    lines.roc(rocobj2, col = "#008600", lwd = 1, ...)
    testobj <- roc.test(rocobj1, rocobj2)
    text(0.5,0.5, labels = paste("p-value =",
     format.pval(testobj$p.value)), adj=c(0, .5))
    legend("bottomright", legend=c(contrasts[coef1], contrasts[coef2]),
     col=c("#1c61b6", "#008600"), lwd=2)
   
    return(list(pval = testobj$p.value,
     AUC1 = auc(rocobj1), pAUC1 = auc(rocobj1,
      partial.auc = c(1, 1 - pAUCthreshold)),
     AUC2 = auc(rocobj2), pAUC1 = auc(rocobj2,
      partial.auc = c(1, 1 - pAUCthreshold)) ))
}