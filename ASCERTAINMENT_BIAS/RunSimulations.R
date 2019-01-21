source("Functions.R")
library(ape)
library(phytools)
library(geiger)

GetGenus <- function(x) {
  return(unname(strsplit(x, "_")[[1]][1]))
}

DoGenusRep <- function(genus, phy, data, ...) {
  taxa <- rownames(data)
  matching.elements <- which(sapply(rownames(data), GetGenus)==genus)
  data <- data[matching.elements]
  names(data) <- taxa[matching.elements]
  pruned <- geiger::treedata(phy, data, sort=TRUE, warnings=FALSE)
  char.corhmm <- data.frame(species=rownames(pruned$data), region=unname(pruned$data), stringsAsFactors=FALSE)
  phy <- pruned$phy
  return(RunAndSummarize(phy, char.corhmm, ...))
}

RunAndSummarize <- function(phy, data, ...) {
  count0 <-  sum(data$region==0)
  count1 <-  sum(data$region==1)
  if(any(c(count0, count1)==0)) {
    rate01 = NA
    rate10 = NA
  } else {
    corhmm.result <- corHMM::corHMM(phy, data, rate.cat=1, node.states="none", root.p="maddfitz", ...)
    rate01 <- corhmm.result$solution[1,2]
    rate10 <- corhmm.result$solution[2,1]
  }
  result.overall <- data.frame(rate01=rate01, rate10=rate10, count0=count0, count1=count1, ntip=Ntip(phy))
  return(result.overall)
}

DoSingleSimRep <- function(phy, rates, mintax=10, maxtax=Inf,...) {
  genera.table <- table(unname(sapply(phy$tip.label, GetGenus)))
  genera.table <- genera.table[which(genera.table>=mintax)]
  genera.table <- genera.table[which(genera.table<=maxtax)]
  genera <- names(genera.table)
  chars <- geiger::sim.char(phy, par=list(rates), nsim=1, model="discrete", root=sample.int(2,1))
  chars <- chars - 1 #start at 0
  results <- sapply(genera, DoGenusRep, phy=phy, data=chars, ...)
  char.corhmm <- data.frame(species=rownames(chars), region=unname(chars), stringsAsFactors=FALSE)
  rownames(char.corhmm) <- char.corhmm$species
  entire.tree <- RunAndSummarize(phy, char.corhmm, ...)
  all.results <- cbind(t(entire.tree),results)
  colnames(all.results)[1] <- "ENTIRE"
  return(all.results)
}


load("ResultThreshold0.5.rda")
phy <- result1_0.5$phy
rates <- result1_0.5$solution
diag(rates) <- -rowSums(rates, na.rm=TRUE)
sim.results.list <- list()
for (i in sequence(100)) {
  sim.results.list[[i]] <- try(DoSingleSimRep(phy, rates, mintax=10, nstarts=4, n.cores=4))
  save(list=ls(), file="SimulationResults.rda")
  try(system("git add SimulationResults.rda"))
  try(system("git commit -m'Sim result' -a"))
}
