library(corHMM)
yearn::yearn(RateViz)


load("EmpiricalResults.rda")
for (i in sequence(length(all.results))) {
  pdf(file=paste0(names(all.results)[[i]], "_Figure_AICc_",all.results[[i]]$AICc,".pdf"), width=10, height=10)
  PlotBubbleMatrix(all.results[[i]]$solution)
  dev.off()
}
