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
  return(RunAndSummarize(phy, char.corhmm, run.type="genus",...))
}

DoCladeRep <- function(clade, rank, phy, data, taxonomy, ...) {
  matching.elements <- which(taxonomy[,which(colnames(taxonomy)==rank)]==clade)
  taxa.to.delete <- taxonomy$gs[-matching.elements]
  phy <- ape::drop.tip(phy, taxa.to.delete)
  data <- data[,1,1]
  pruned <- geiger::treedata(phy, data, sort=TRUE, warnings=FALSE)
  char.corhmm <- data.frame(species=rownames(pruned$data), region=unname(pruned$data), stringsAsFactors=FALSE)
  phy <- pruned$phy
  return(RunAndSummarize(phy, char.corhmm, run.type=rank,...))
}

DoRandomRep <- function(ntax, phy, data, run.type=NULL, ...) {
  phy <- geiger::drop.random(phy, ape::Ntip(phy)-ntax)
  data <- data[,1,1]
  pruned <- geiger::treedata(phy, data, sort=TRUE, warnings=FALSE)
  char.corhmm <- data.frame(species=rownames(pruned$data), region=unname(pruned$data), stringsAsFactors=FALSE)
  phy <- pruned$phy
  if(is.null(run.type)) {
    run.type <- paste0("random_",ntax,"_taxa")
  }
  return(RunAndSummarize(phy, char.corhmm, run.type=run.type,...))
}

RunAndSummarize <- function(phy, data, run.type, ...) {
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
  result.overall <- data.frame(rate01=rate01, rate10=rate10, count0=count0, count1=count1, ntip=Ntip(phy), depth=max(ape::branching.times(phy)), total.length = sum(phy$edge.length), label=run.type, stringsAsFactors=FALSE)
  return(result.overall)
}

DoSingleSimRep <- function(phy, rates, mintax=50, maxtax=5000, taxonomy,...) {
  genera.table <- table(unname(sapply(phy$tip.label, GetGenus)))
  genera.table <- genera.table[which(genera.table>=mintax)]
  genera.table <- genera.table[which(genera.table<=maxtax)]
  genera <- names(genera.table)

  families.table <- table(unname(taxonomy$family))
  families.table <- families.table[which(families.table>=mintax)]
  families.table <- families.table[which(families.table<=maxtax)]
  families <- names(families.table)

  orders.table <- table(unname(taxonomy$order))
  orders.table <- orders.table[which(orders.table>=mintax)]
  orders.table <- orders.table[which(orders.table<=maxtax)]
  orders <- names(orders.table)

  chars <- geiger::sim.char(phy, par=list(rates), nsim=1, model="discrete", root=2) #2 = tropical when they're on a 1,2 scale
  chars <- chars - 1 #go back to 0,1 scale
  results.genus <- sapply(genera, DoGenusRep, phy=phy, data=chars, ...)
  results.family <- sapply(families, DoCladeRep, phy=phy, data=chars, taxonomy=taxonomy, rank="family",...)
  results.order <- sapply(orders, DoCladeRep, phy=phy, data=chars, taxonomy=taxonomy, rank="order",...)



  ntax.vector <- c(unlist(results.genus["ntip",]), unlist(results.family["ntip",]), unlist(results.order["ntip",]))
  names(ntax.vector) <- paste0(names(ntax.vector), "_RandomTwin")

  results.random <- sapply(ntax.vector, DoRandomRep, phy=phy, data=chars, ...)

#  char.corhmm <- data.frame(species=rownames(chars), region=unname(chars), stringsAsFactors=FALSE)
#rownames(char.corhmm) <- char.corhmm$species
  #entire.tree <- RunAndSummarize(phy, char.corhmm, run.type="entire", ...)
  all.results <- cbind(results.genus, results.family, results.order, results.random)
  #colnames(all.results)[1] <- "ENTIRE"
  return(all.results)
}

taxonomy <- read.csv("Campanulid.reference.table.csv", stringsAsFactors=FALSE)
taxonomy <- taxonomy[which(nchar(taxonomy$gs)>0),]

load("ResultThreshold0.5.rda")
phy <- result1_0.5$phy
rates <- result1_0.5$solution
diag(rates) <- -rowSums(rates, na.rm=TRUE)
sim.results.list <- list()
for (i in sequence(100)) {
  sim.results.list[[i]] <- try(DoSingleSimRep(phy, rates, taxonomy=taxonomy, mintax=50, maxtax=5000, nstarts=3, n.cores=3))
  save(list=ls(), file="SimulationResultsPaired.rda")
  try(system("git add SimulationResultsPaired.rda"))
  try(system("git commit -m'Sim result' -a"))
  try(system("git push"))
}
