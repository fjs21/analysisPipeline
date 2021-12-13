
ROC.compare2gs <- 
function (gs1, gs2, fit3 = fit3, coef, main, pAUCthreshold = 0.2, add = FALSE, 
    gs1.name = NULL, gs2.name = NULL, ...) 
{
    require(verification)
    if (class(gs1) == "GeneSet") {
        require(GSEABase)
        gs1.EG <- geneIds(gs1)
        gs1.name <- setName(gs1)		
    }
    else {
        gs1.EG <- as.character(gs1)
        if(is.null(gs1.name)) {gs1.name <- "Set 1"}
    }
    if (class(gs2) == "GeneSet") {
        require(GSEABase)
        gs2.EG <- geneIds(gs2)
    }
    else {
        gs2.EG <- as.character(gs2)
        if(is.null(gs2.name)) {gs2.name <- "Set 2"}
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

    truth1 <- (EG %in% gs1.EG) * 1
    require(pROC)
    rocobj1 <- roc(truth1, data, direction = "<")
    plot.roc(rocobj1, col = "#1c61b6", identity.col = 'red',
     lwd = 1, main = main, ...)
    abline(v = 1 - pAUCthreshold, col = "black", lty = 2)

    truth2 <- (EG %in% gs2.EG) * 1
    rocobj2 <- roc(truth2, data, direction = "<")
    lines.roc(rocobj2, col = "#008600", lwd = 1, ...)

    testobj <- roc.test(rocobj1, rocobj2)
    text(0.5,0.5, labels = paste("p-value =",
     format.pval(testobj$p.value)), adj=c(0, .5))
    legend("bottomright", legend=c(gs1.name, gs2.name),
     col=c("#1c61b6", "#008600"), lwd=2)
   
    return(list(pval = testobj$p.value,
     AUC1 = auc(rocobj1), pAUC1 = auc(rocobj1,
      partial.auc = c(1, 1 - pAUCthreshold)),
     AUC2 = auc(rocobj2), pAUC2 = auc(rocobj2,
      partial.auc = c(1, 1 - pAUCthreshold)) ))
}