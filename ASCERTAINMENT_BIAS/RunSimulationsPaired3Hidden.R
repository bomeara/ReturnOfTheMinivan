rm(list=ls())
source("Functions.R")
library(ape)
library(phytools)
library(geiger)
load("EmpiricalResults.rda")

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
    #data <- data[,1,1] #Don't need this anymore, as do earlier
    pruned <- geiger::treedata(phy, data, sort=TRUE, warnings=FALSE)
    char.corhmm <- data.frame(species=rownames(pruned$data), region=unname(pruned$data), stringsAsFactors=FALSE)
    phy <- pruned$phy
    return(RunAndSummarize(phy, char.corhmm, run.type=rank,...))
}

DoRandomRep <- function(ntax, phy, data, run.type=NULL, ...) {
    phy <- geiger::drop.random(phy, ape::Ntip(phy)-ntax)
    #data <- data[,1,1] #Don't need this anymore, as do earlier
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

DoSingleSimRep <- function(phy, rates, mintax=50, maxtax=5000, taxonomy, count.index, ...) {
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
    
    chars.original <- geiger::sim.char(phy, par=list(rates), nsim=1, model="discrete", root=sample(c(2,4,6),1))[,1,1] #2 = tropical when they're on a 1,2 scale
    chars.transformed <- matrix((chars.original-1)%%2, ncol=1) #to get back to zero and 1
    rownames(chars.transformed) <- names(chars.original)
    results.genus <- sapply(genera, DoGenusRep, phy=phy, data=chars.transformed, ...)
    results.family <- sapply(families, DoCladeRep, phy=phy, data=chars.transformed, taxonomy=taxonomy, rank="family",...)
    results.order <- sapply(orders, DoCladeRep, phy=phy, data=chars.transformed, taxonomy=taxonomy, rank="order",...)
    
    
    
    ntax.vector <- c(unlist(results.genus["ntip",]), unlist(results.family["ntip",]), unlist(results.order["ntip",]))
    names(ntax.vector) <- paste0(names(ntax.vector), "_RandomTwin")
    
    results.random <- sapply(ntax.vector, DoRandomRep, phy=phy, data=chars.transformed, ...)
    
    #  char.corhmm <- data.frame(species=rownames(chars), region=unname(chars), stringsAsFactors=FALSE)
    #rownames(char.corhmm) <- char.corhmm$species
    #entire.tree <- RunAndSummarize(phy, char.corhmm, run.type="entire", ...)
    all.results <- cbind(results.genus, results.family, results.order, results.random)
    #colnames(all.results)[1] <- "ENTIRE"
    save(chars.original, chars.transformed, phy, all.results, file=paste0("IndividualSimulationResultsPaired3Hidden_replicate", count.index, ".rda"))
    
    return(all.results)
}

taxonomy <- read.csv("Campanulid.reference.table.csv", stringsAsFactors=FALSE)
taxonomy <- taxonomy[which(nchar(taxonomy$gs)>0),]


rates <- all.results[["RateCat_3_Threshold_0.5"]]$solution
rates[which(is.na(rates))] <- 0
diag(rates) <- -rowSums(rates, na.rm=TRUE)
sim.results.list <- list()
for (i in sequence(100)) {
    sim.results.list[[i]] <- try(DoSingleSimRep(phy, rates, taxonomy=taxonomy, count.index=i, mintax=50, maxtax=5000, nstarts=3, n.cores=3))
    save(list=ls(), file="SimulationResultsPaired3Hidden.rda")
    try(system("git add SimulationResultsPaired3Hidden.rda"))
    try(system("git commit -m'Sim result' -a"))
    try(system("git push"))
}
