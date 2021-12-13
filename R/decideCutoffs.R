decideCutoffs <-
function(fit3, fc = 2) {

 require(limma)
 
 decisionM <- matrix(c("fdr", 0.01, "fdr", 0.05, "fdr", 0.1,
  "none", 0.001, "none", 0.005, "none", 0.01),nrow = 6 ,ncol = 2,byrow = TRUE)

 testSummarys = as.list(1:nrow(decisionM))
 testSigNumber = 0
 maxmax = 0

 for (test in 1:nrow(decisionM)) {
   testSummarys[[test]] = list()
   results <- decideTests(fit3, adjust.method = decisionM[test,1], 
    p.value = as.numeric(decisionM[test,2]))
   # summarize results
   testSummarys[[test]]$all <- summary(results)
   # count number of features with a significantly UP regulated result
   # why slow?
   testSigNumber[test] <- length(which(apply(results,1,function(x) any(x != 0,na.rm = TRUE))))
   # find maximum value of regulated features
   maxmax <- max(c(maxmax, as.vector(testSummarys[[test]]$all[c(1,3),])))
   results <- decideTests(fit3, adjust.method = decisionM[test,1], 
    p.value = as.numeric(decisionM[test,2]), lfc = log2(fc) )
   testSummarys[[test]]$fc <- summary(results)
 }

 par(mfrow=c(2,nrow(decisionM)/2))
 for (test in 1:nrow(decisionM)) {
  UP = -testSummarys[[test]]$all[3,]
  DOWN = +testSummarys[[test]]$all[1,]
  
  barplot2(UP, horiz = TRUE, col="red", xlim = c(-maxmax,maxmax), names.arg = NULL,
   main = paste("Summary - p<", decisionM[test,2], ", ", decisionM[test,1],
   " (", testSigNumber[test], ")", sep = "")) 
  barplot2(DOWN, horiz = TRUE, col="green", names.arg = NULL, add = TRUE) -> yy

  xx <- vector("integer",ncol(testSummarys[[test]]$all))
  text(xx,yy, colnames(testSummarys[[test]]$all))
  text(min(UP)  , yy-0.1, format(abs(UP),digits=3), adj=c(1,0.5), xpd = TRUE)
  text(max(DOWN), yy+0.1, format(abs(DOWN),digits=3), adj=c(0,0.5), xpd = TRUE )

  UP = -testSummarys[[test]]$fc[3,]
  DOWN = +testSummarys[[test]]$fc[1,]
  
  barplot2(UP, horiz = TRUE, col="orange", names.arg = NULL, add=TRUE) 
  barplot2(DOWN, horiz = TRUE, col="lightgreen", names.arg = NULL, add=TRUE) -> yy

  xx <- vector("integer",ncol(testSummarys[[test]]$fc))
  text(xx,yy, colnames(testSummarys[[test]]$fc))
  text(UP  , rep(yy+0.1,length(UP)), format(abs(UP),digits=3), adj=c(1,0.5) )
  text(DOWN, rep(yy-0.1,length(DOWN)), format(abs(DOWN),digits=3), adj=c(0,0.5) )  
 }

 legend("bottomright",
  c("UP","DOWN",paste("UP >", fc, " fold", sep=""),
   paste("DOWN >", fc, " fold", sep="")), cex = 0.7,
  fill = c("red","green","orange","lightgreen"), bg = "white", ncol = 2 )

}

