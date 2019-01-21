library(corHMM)
load("SimulationResultsMultipleLevels.rda")
library(viridis)
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

all.sims <- sapply(sequence(ncol(sim.results.list[[1]])), SummarizeIndividualSim, sim.results.list=sim.results.list)

all.sims.df <- data.frame(all.sims[,1], stringsAsFactors=FALSE)
for (i in 2:ncol(sim.results.list[[1]])) {
  all.sims.df <- rbind(all.sims.df, data.frame(all.sims[,i], stringsAsFactors=FALSE))
}

pdf(file="Rates.pdf", width=10, height=10)
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


pdf(file="AverageRateVsNtax.pdf", width=10, height=10)
plot(x=range(all.sims.df$ntip.median), y=range(all.sims.df$average.rate.median), bty="n", xlab="Number of taxa", ylab="Average transition rate", type="n", log="xy")
abline(h=(truerate01+truerate10)/2, col="red", lty="dotted")
points(all.sims.df.taxa$ntip.median, all.sims.df.taxa$average.rate.median, pch=19)
points(all.sims.df.random$ntip.median, all.sims.df.random$average.rate.median, pch=19, col="red")
dev.off()


pdf(file="RateRatioVsNtax.pdf", width=10, height=10)
plot(x=range(all.sims.df$ntip.median), y=range(all.sims.df$tropic.to.temperate.fraction.median), bty="n", xlab="Number of taxa", ylab="Tropical to temperate rate fraction", type="n", log="x")
abline(h=truerate10/(truerate01+truerate10), col="red", lty="dotted")
points(all.sims.df.taxa$ntip.median, all.sims.df.taxa$tropic.to.temperate.fraction.median, pch=19)
points(all.sims.df.random$ntip.median, all.sims.df.random$tropic.to.temperate.fraction.median, pch=19, col="red")
dev.off()


write.csv(all.sims.df, file="AllSimRates.csv")
