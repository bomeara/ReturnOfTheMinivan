source("Functions.R")
tree <- ape::read.tree("camp.fossils.rootcon.dated.2.tre")
print(tree)
data <- read.delim("camp.trop.temp.txt", stringsAsFactors = FALSE)
head(data)
rownames(data) <- data[,1]
pruned <- geiger::treedata(tree, data, warnings=TRUE)
tree <- pruned$phy
data <- pruned$data
data <- data.frame(species=data[,"species"], region=as.numeric(data[,"region"]), tot.records=as.numeric(data[,"tot.records"]), stringsAsFactors=FALSE)
head(data)
print(tree)
rate.cats <- c(1,2,3)
thresholds <- c(0.5, 0.05, 0.95)
all.results <- list()
for (rate.cat.index in sequence(length(rate.cats))) {
  for(threshold.index in sequence(length(thresholds))) {
    print(paste0("Starting RateCat_",rate.cats[rate.cat.index], "_Threshold_",thresholds[threshold.index]))
    local.result <- DoSingleAnalysis(tree, data, thresholds[threshold.index], n.cores=2, rate.cat=rate.cats[rate.cat.index], nstarts=2, node.states="none")
    all.results[[length(all.results)+1]] <- local.result
    names(all.results)[length(all.results)] <- paste0("RateCat_",rate.cats[rate.cat.index], "_Threshold_",thresholds[threshold.index])
    save(list=ls(), file="EmpiricalResults.rda")
  }
}
