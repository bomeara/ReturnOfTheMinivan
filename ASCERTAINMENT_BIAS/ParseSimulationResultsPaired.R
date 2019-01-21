rm(list=ls())
library(corHMM)
load("SimulationResultsPaired.rda")
library(viridis)
library(car)
truerate01 <- rates[1,2]
truerate10 <- rates[2,1]

#remember 0 to 1 is temperate to tropical

# First, need to label the random taxa
for (i in sequence(length(sim.results.list))) {
  local.data <- sim.results.list[[i]]
  missing.labels <- which(nchar(colnames(local.data))==0)
  colnames(local.data)[missing.labels] <- local.data$label[missing.labels]
  sim.results.list[[i]] <- local.data
}

SummarizeIndividualSim <- function(position, sim.results.list) {
  sim.kind <- sim.results.list[[1]]["label", position][[1]]
  sim.name <- colnames(sim.results.list[[1]])[position]
  sub.list <- list()
  for(i in sequence(length(sim.results.list))) {
    sub.list[[i]] <- sim.results.list[[i]][,position]
    sub.list[[i]] <- sub.list[[i]][-which(names(sub.list[[i]])=="label")]
  }
  final.df <- t(simplify2array(sub.list))
  final.df <- data.frame(apply(final.df, 2, as.numeric))
  final.df$average.rate <- (final.df$rate01 + final.df$rate10)/2
  final.df$tropic.to.temperate.fraction = final.df$rate10 / (final.df$rate01 + final.df$rate10)
  medians <- apply(final.df, 2, median, na.rm=TRUE)
  sds <- apply(final.df, 2, sd, na.rm=TRUE)
  names(medians) <- paste0(names(medians), ".median")
  names(sds) <- paste0(names(sds), ".sd")
  invariant <- data.frame(invariant.proportion = sum(is.na(final.df$rate01)))
  result <- cbind(data.frame(kind=sim.kind, name=sim.name, stringsAsFactors=FALSE), data.frame(t(medians)), data.frame(t(sds)), invariant)
  return(result)
}

ComputeRateMeasures <- function(m, tr01=truerate01, tr10=truerate10) {
    for(i in sequence(ncol(m)-1)) { #last col is label
        m[,i] <- as.numeric(m[,i])
    }

    m$tropic.to.temperate.fraction <- m$rate10 / (m$rate01 + m$rate10)
    m$average.rate <- (m$rate01 + m$rate10)/2
    m$tropic.to.temperate.fraction.deviation <- m$tropic.to.temperate.fraction - tr10/(tr01+tr10)
    m$abs.tropic.to.temperate.fraction.deviation <- abs(m$tropic.to.temperate.fraction.deviation)
    m$average.rate.deviation <- m$average.rate - (tr01 + tr01)/2
    m$abs.average.rate.deviation <- abs(m$average.rate.deviation)
    return(m)
}

GetPairsFromIndividualRuns <- function(m) {
  empirical <- ComputeRateMeasures(data.frame(t(m[,which(!grepl("RandomTwin", colnames(m)))]), stringsAsFactors = FALSE))
  empirical$clade.name <- rownames(empirical)
  colnames(empirical) <- paste0("empirical.", colnames(empirical))
  random <- ComputeRateMeasures(data.frame(t(m[,which(grepl("RandomTwin", colnames(m)))]), stringsAsFactors = FALSE))
  rownames(random) <- gsub("_RandomTwin", "", rownames(random))
  random$clade.name <- rownames(random)
  colnames(random) <- paste0("random.", colnames(random))
  final.df <- merge(empirical, random, by="row.names", all.x=TRUE)
  rownames(final.df) <- final.df$Row.names
  final.df <- final.df[,-1]
  return(final.df)
}

GetModelResults <- function(model) {
    local.summary <- summary(model)
    local.conf <- confint(model)
    result.df <- data.frame(intercept=local.summary$coefficients[1,1], slope=local.summary$coefficients[2,1], intercept.p=local.summary$coefficients[1,4], slope.p=local.summary$coefficients[2,4], intercept.025=local.conf[1,1], intercept.975=local.conf[1,2], slope.025=local.conf[2,1], slope.975=local.conf[2,2], r.squared=local.summary$r.squared, adj.r.squared=local.summary$adj.r.squared)
    return(result.df)
}

RunModels <- function(pairs.df) {

    formulas <- c(
        empirical.average.rate ~ log(empirical.ntip),
        random.average.rate ~ log(random.ntip),
        empirical.tropic.to.temperate.fraction ~ log(empirical.ntip),
        random.tropic.to.temperate.fraction ~ log(random.ntip),
        abs(empirical.average.rate.deviation) ~ log(empirical.ntip),
        abs(random.average.rate.deviation) ~ log(random.ntip),
        abs(empirical.tropic.to.temperate.fraction.deviation) ~ log(empirical.ntip),
        abs(random.tropic.to.temperate.fraction.deviation) ~ log(random.ntip)
    )
    names(formulas) <- c(
        "Empirical.rate",
        "Random.rate",
        "Empirical.tropic",
        "Random.tropic",
        "Empirical.rate.deviation",
        "Random.rate.deviation",
        "Empirical.tropic.deviation",
        "Random.tropic.deviation"
    )

    model.returns <- lapply(formulas, lm, data=pairs.df)
    model.summaries <- lapply(model.returns, GetModelResults)
    model.summaries <- do.call(rbind, model.summaries)
    # #plot(pairs.df$empirical.ntip, pairs.df$empirical.average.rate, pch=19, col=rgb(0,0,0,.5), log="xy")
    # model <- lm(empirical.average.rate ~ log(empirical.ntip), data=pairs.df)
    # print(model)
    # print(summary(model))
    # print(confint(model))
    #
    # #plot(pairs.df$empirical.ntip, pairs.df$empirical.average.rate.deviation, pch=19, col=rgb(0,0,0,.5), log="xy")
    # model <- lm(abs(empirical.average.rate.deviation) ~ log(empirical.ntip), data=pairs.df)
    # print(model)
    # print(summary(model))
    # print(confint(model))
    #
    # model <- lm(abs(random.average.rate.deviation) ~ log(random.ntip), data=pairs.df)
    # print(model)
    # print(summary(model))
    # print(confint(model))
    #
    # #plot(pairs.df$empirical.ntip, pairs.df$empirical.tropic.to.temperate.fraction.deviation, pch=19, col=rgb(0,0,0,.5), log="xy")
    # model <- lm(abs(empirical.tropic.to.temperate.fraction.deviation) ~ log(empirical.ntip), data=pairs.df)
    # print(model)
    # print(summary(model))
    # print(confint(model))
    #
    # model <- lm(abs(random.tropic.to.temperate.fraction.deviation) ~ log(random.ntip), data=pairs.df)
    # print(model)
    # print(summary(model))
    # print(confint(model))
}

DoPairedTest <- function(trait, pairs.df) {
    x1 <- pairs.df[, which(colnames(pairs.df)==paste0("empirical.", trait))]
    x2 <- pairs.df[, which(colnames(pairs.df)==paste0("random.", trait))]
    w.result <- wilcox.test(x1, x2, paired=TRUE)
    t.result <- t.test(x1, x2, paired=TRUE)
    fligner.data <- data.frame(x=c(x1, x2), y=c(rep(0, length(x1)), rep(1, length(x2))))
    fligner.result <- stats::fligner.test(x~y, data=fligner.data)
    return(data.frame(wilcox.p = w.result$p.value, fligner.p = fligner.result$p.value, ttest.p = t.result$p.value, ttest.025=t.result$conf.int[1], ttest.975=t.result$conf.int[2]))
}

DoAllPairedTests <- function(pairs.df) {
     traits <- c("average.rate", "tropic.to.temperate.fraction", "abs.average.rate.deviation", "abs.tropic.to.temperate.fraction.deviation")
     return(data.frame(t(sapply(traits, DoPairedTest, pairs.df=pairs.df))))
}

PlotPairs <- function(pairs.df, trait, trait.name, main, ...) {
    pairs.df <- pairs.df[order(runif(nrow(pairs.df))),] #reorder randomly
    x1 <- pairs.df[, which(colnames(pairs.df)==paste0("empirical.", trait))]
    x2 <- pairs.df[, which(colnames(pairs.df)==paste0("random.", trait))]
    plot(x=c(-.25,1.25), y=range(c(x1,x2)), bty="n", xaxt="n", xlab="", ylab=trait.name, main=paste0(main,'\nMedians = ',signif(median(x1),2), ' and ', signif(median(x2),2), '\nSDs = ', signif(sd(x1),2), ' and ', signif(sd(x2),2)), cex.main=0.8, type="n", ...)
    text(0, median(x1), "Named Clades", pos=2, cex=0.8, srt=90, adj=c(0.5, 0.5))
    text(1, median(x2), "Random Subset", pos=4, cex=0.8, srt=-90)
    mtext("(A)",side=3, line=0, adj=0)
    for(i in sequence(length(x1))) {
        color <- rgb(1,0,0,.1)
        if(x1[i]<x2[i]) {
            color <- rgb(0,0,1,.1)
        }
        lines(c(0,1), c(x1[i], x2[i]), col=color)
    }
    points(0, median(x1), pch=19, col="black")
    points(1, median(x2), pch=19, col="black")
    lines(c(0,1), c(median(x1), median(x2)), col="black")
    #arrows(x0=0, y0=max(min(x1),median(x1) - sd(x1)), x1=0, y1=median(x1) + sd(x1), angle=90, col="black", length=0.1, lwd=2, code=3) #min to deal with log scale issues
    #arrows(x0=1, y0=max(min(x2),median(x2) - sd(x2)), x1=1, y1=median(x2) + sd(x2), angle=90, col="black", length=0.1, lwd=2, code=3) #min to deal with log scale issues

}

#Do all sims individually


list.of.pairs <- lapply(sim.results.list, GetPairsFromIndividualRuns)
for (i in sequence(length(list.of.pairs))) {
    list.of.pairs[[i]]$replicate <- i
}

pairs.df <- do.call(rbind, list.of.pairs)

pairs.df <- subset(pairs.df, pairs.df$empirical.ntip>=50) #some taxa had too few observations

pairs.df <- subset(pairs.df, !is.na(pairs.df$empirical.rate01)) # no NA

print("All sims where there were at least some variation")





at.least.3.taxa.each.state <- subset(pairs.df, pairs.df$empirical.count0>=3)
at.least.3.taxa.each.state <- subset(at.least.3.taxa.each.state, at.least.3.taxa.each.state$empirical.count1>=3)
at.least.3.taxa.each.state <- subset(at.least.3.taxa.each.state, at.least.3.taxa.each.state$random.count0>=3)
at.least.3.taxa.each.state <- subset(at.least.3.taxa.each.state, at.least.3.taxa.each.state$random.count1>=3)

print("All sims where there there were at least 10% of taxa in each state")
at.least.10.percent.each.state <- subset(pairs.df, pairs.df$empirical.count0/pairs.df$empirical.ntip >= 0.1)
at.least.10.percent.each.state <- subset(at.least.10.percent.each.state, at.least.10.percent.each.state$empirical.count1/at.least.10.percent.each.state$empirical.ntip >= 0.1)
at.least.10.percent.each.state <- subset(at.least.10.percent.each.state, at.least.10.percent.each.state$random.count1/at.least.10.percent.each.state$random.ntip >= 0.1)
at.least.10.percent.each.state <- subset(at.least.10.percent.each.state, at.least.10.percent.each.state$random.count1/at.least.10.percent.each.state$random.ntip >= 0.1)


datasets <- list(pairs.df, at.least.3.taxa.each.state, at.least.10.percent.each.state)
names(datasets) <- c("variation.some", "at.least.3.taxa.each.state", "at.least.10.percent.taxa.each.state")

datasets <- c(datasets, lapply(datasets, subset, empirical.ntip>=200))
names(datasets)[4:6] <- paste0(names(datasets)[1:3], ".200plus.taxa.total")

correlations <- lapply(datasets, RunModels)
paired.tests <- lapply(datasets, DoAllPairedTests)

for (i in sequence(length(correlations))) {
    correlations[[i]]$dataset <- names(correlations)[i]
    correlations[[i]]$npoints <- nrow(datasets[[i]])
    paired.tests[[i]]$dataset <- names(paired.tests)[i]
    paired.tests[[i]]$npoints <- nrow(datasets[[i]])
}

save(correlations, paired.tests, file="OrganizedSummaryOfPairedSims.rda")

trait.names.pretty <- c("Average transition rate", "Tropical to temperate rate proportion", "Error in transition rate", "Error in tropical rate proportion")
logs <- c("y", "", "y", "y")

for (dataset.index in sequence(length(datasets))) {
    pdf(file=paste0("Figure_", names(datasets)[dataset.index], "_pairedtests.pdf"), width=10, height=10)
    par(mfcol=c(2,2))
    #par(mar=c(9,4,4,2)+.1)
    for (trait.index in sequence(nrow(paired.tests[[dataset.index]]))) {
      #PlotPairs <- function(pairs.df, trait, trait.name, main, ...)
        PlotPairs(datasets[[dataset.index]], trait=rownames(paired.tests[[dataset.index]])[trait.index], trait.name=trait.names.pretty[trait.index], main=paste0("Wilcoxon p = ",  signif(as.numeric((paired.tests[[dataset.index]])$wilcox.p[trait.index]),2), '\nFligner p = ',  signif(as.numeric((paired.tests[[dataset.index]])$fligner.p[trait.index]),2)), log=logs[trait.index])
    }
    dev.off()
    Sys.sleep(2)
    #system(paste0("open Figure_", names(datasets)[dataset.index], "_pairedtests.pdf"))
}

all.clade.names <- unique(pairs.df$empirical.clade.name)
for (clade.index in sequence(length(all.clade.names))) {
  pdf(file=paste0("Figure_", all.clade.names[clade.index], "_pairedtests.pdf"), width=7.5, height=7.5)
  par(mfcol=c(2,2))
  pairs.taxon <- subset(pairs.df, empirical.clade.name==all.clade.names[clade.index])
  paired.tests.taxon <- DoAllPairedTests(pairs.taxon)
  #par(mar=c(9,4,4,2)+.1)
  for (trait.index in sequence(nrow(paired.tests.taxon))) {
    #PlotPairs <- function(pairs.df, trait, trait.name, main, ...)
    PlotPairs(pairs.taxon, trait=rownames(paired.tests.taxon)[trait.index], trait.name=trait.names.pretty[trait.index], main=paste0(all.clade.names[clade.index],'\n(',pairs.taxon$empirical.ntip[1], " taxa sampled)","\nWilcoxon p = ",  signif(as.numeric(paired.tests.taxon$wilcox.p[trait.index]),2), '\ntFligner p = ',  signif(as.numeric(paired.tests.taxon$fligner.p[trait.index]),2)), log=logs[trait.index], cex.main=0.8)
  }
  dev.off()
  Sys.sleep(2)
  #system(paste0("open Figure_", all.clade.names[clade.index], "_pairedtests.pdf"))
}










# Below is summarizing across sims.

all.sims <- sapply(sequence(ncol(sim.results.list[[1]])), SummarizeIndividualSim, sim.results.list=sim.results.list)

all.sims.df <- data.frame(all.sims[,1], stringsAsFactors=FALSE)
for (i in 2:ncol(sim.results.list[[1]])) {
  all.sims.df <- rbind(all.sims.df, data.frame(all.sims[,i], stringsAsFactors=FALSE))
}

all.sims.df <- all.sims.df[order(all.sims.df$name),]

all.sims.df <- all.sims.df[which(all.sims.df$ntip.median>=50),]


pdf(file="RatesPaired.pdf", width=10, height=10)
plot(x=range(all.sims.df$average.rate.median), y=range(all.sims.df$tropic.to.temperate.fraction.median), bty="n", type="n", xlab="Average transition rate", ylab="Tropical to temperate rate fraction", log="x")

color.fn <- colorRamp(viridis(256))

abline(h=truerate10/(truerate01+truerate10), col="red", lty="dotted")
abline(v=(truerate01+truerate10)/2, col="red", lty="dotted")

arrows(x0=all.sims.df$average.rate.median, y0=all.sims.df$tropic.to.temperate.fraction.median + all.sims.df$tropic.to.temperate.fraction.sd, x1=all.sims.df$average.rate.median, y1=all.sims.df$tropic.to.temperate.fraction.median - all.sims.df$tropic.to.temperate.fraction.sd, length=0, col="gray")

arrows(x0=all.sims.df$average.rate.median - all.sims.df$average.rate.sd, y0=all.sims.df$tropic.to.temperate.fraction.median, x1=all.sims.df$average.rate.median + all.sims.df$average.rate.sd, y1=all.sims.df$tropic.to.temperate.fraction.median, length=0, col="gray")


points(all.sims.df$average.rate.median, all.sims.df$tropic.to.temperate.fraction.median, pch=19, col=rgb(color.fn(sqrt(all.sims.df$ntip.median/max(all.sims.df$ntip.median))), maxColorValue=256))
dev.off()

all.sims.df.genus <- subset(all.sims.df, kind=="genus")
all.sims.df.order <- subset(all.sims.df, kind=="order")
all.sims.df.family <- subset(all.sims.df, kind=="family")
all.sims.df.taxa <- rbind(all.sims.df.genus, all.sims.df.order, all.sims.df.family)
all.sims.df.random <- all.sims.df[grepl("random", all.sims.df$kind),]

all.sims.df.paired <- data.frame(taxa=all.sims.df.taxa$name, ntip=all.sims.df.taxa$ntip.median, empirical.average.rate = all.sims.df.taxa$average.rate.median, empirical.tropic.to.temperate.fraction.median = all.sims.df.taxa$tropic.to.temperate.fraction.median,  random.average.rate = NA, random.tropic.to.temperate.fraction.median = NA, stringsAsFactors = FALSE)
for (i in sequence(length(all.sims.df.paired$taxa))) {
    matching.line <- which(paste0(all.sims.df.paired$taxa[i], "_RandomTwin")==all.sims.df.random$name)
    all.sims.df.paired$random.average.rate[i] <- all.sims.df.random$average.rate.median[matching.line]
    all.sims.df.paired$random.tropic.to.temperate.fraction.median[i] <- all.sims.df.random$tropic.to.temperate.fraction.median[matching.line]
}

all.sims.df.paired$empirical.average.rate.deviation <- all.sims.df.paired$empirical.average.rate - (truerate01+truerate10)/2
all.sims.df.paired$random.average.rate.deviation <- all.sims.df.paired$random.average.rate - (truerate01+truerate10)/2
all.sims.df.paired$empirical.tropic.to.temperate.fraction.median.deviation <- all.sims.df.paired$empirical.tropic.to.temperate.fraction.median - truerate10/(truerate01+truerate10)
all.sims.df.paired$random.tropic.to.temperate.fraction.median.deviation <- all.sims.df.paired$random.tropic.to.temperate.fraction.median - truerate10/(truerate01+truerate10)

data.for.wilcox <- data.frame(average.rate.deviation=abs(c(all.sims.df.paired$random.average.rate.deviation, all.sims.df.paired$empirical.average.rate.deviation)), tropic.to.temperate.fraction.median.deviation=abs(c(all.sims.df.paired$random.tropic.to.temperate.fraction.median.deviation, all.sims.df.paired$empirical.tropic.to.temperate.fraction.median.deviation)), random1 = c(rep(1, nrow(all.sims.df.paired)), rep(0, nrow(all.sims.df.paired))))


print("Average rate")
print(t.test(abs(all.sims.df.paired$empirical.average.rate.deviation), abs(all.sims.df.paired$random.average.rate.deviation), paired=TRUE))
print(wilcox.test(average.rate.deviation ~ random1, data=data.for.wilcox))

print("Tropical temperate rate")
print(t.test(abs(all.sims.df.paired$random.average.rate.deviation), abs(all.sims.df.paired$random.tropic.to.temperate.fraction.median.deviation), paired=TRUE))
print(wilcox.test(tropic.to.temperate.fraction.median.deviation ~ random1, data=data.for.wilcox))


pdf(file="AverageRateVsNtaxPaired.pdf", width=10, height=10)
plot(x=range(all.sims.df$ntip.median), y=range(all.sims.df$average.rate.median), bty="n", xlab="Number of taxa", ylab="Average transition rate", type="n", log="xy")
abline(h=(truerate01+truerate10)/2, col="red", lty="dotted")
points(all.sims.df.taxa$ntip.median, all.sims.df.taxa$average.rate.median, pch=19)
points(all.sims.df.random$ntip.median, all.sims.df.random$average.rate.median, pch=19, col="red")
for (i in sequence(nrow(all.sims.df.paired))) {
    lines(x=rep(all.sims.df.paired$ntip[i],2), y=c(all.sims.df.paired$empirical.average.rate[i], all.sims.df.paired$random.average.rate[i]))
}
dev.off()


pdf(file="RateRatioVsNtaxPaired.pdf", width=10, height=10)
plot(x=range(all.sims.df$ntip.median), y=range(all.sims.df$tropic.to.temperate.fraction.median), bty="n", xlab="Number of taxa", ylab="Tropical to temperate rate fraction", type="n", log="x")
abline(h=truerate10/(truerate01+truerate10), col="red", lty="dotted")
points(all.sims.df.taxa$ntip.median, all.sims.df.taxa$tropic.to.temperate.fraction.median, pch=19)
points(all.sims.df.random$ntip.median, all.sims.df.random$tropic.to.temperate.fraction.median, pch=19, col="red")
dev.off()


write.csv(all.sims.df, file="AllSimRates.csv")
