library(data.table)

SummarizeBiomes <- function(ecoregions.set){
    res <- c()
    names(ecoregions.set) <- c("species", "biome", "lat")
    species.names <- unique(ecoregions.set$species)
    for(species.index in 1:length(species.names)){
        print(species.index)
        subset.mat <- ecoregions.set[which(ecoregions.set$species==species.names[species.index]),]
        #Step 1: Remove the -9999's, 98's, 99's.
        subset.mat <- subset.mat[which(!subset.mat$biome==-9999),]
        subset.mat <- subset.mat[which(!subset.mat$biome==98),]
        subset.mat <- subset.mat[which(!subset.mat$biome==99),]
        #Step 2: Define tropical vs. temperate:
        tmp.vector <- numeric(dim(subset.mat)[1])
        if(dim(subset.mat)[1]>0){
            for(record.index in 1:dim(subset.mat)[1]){
                if(subset.mat$biome[record.index] == 1 | subset.mat$biome[record.index] == 2 | subset.mat$biome[record.index] == 3 |subset.mat$biome[record.index] == 7 | subset.mat$biome[record.index] == 14){
                    tmp.vector[record.index] <- 1
                }
                if(subset.mat$biome[record.index] == 10){
                    if(subset.mat$lat[record.index] < 23.5 & subset.mat$lat[record.index] > -23.5){
                        tmp.vector[record.index] <- 1
                    }
                }
            }
        }
        if(dim(subset.mat)[1]>0){
            final.tally <- sum(tmp.vector)/ length(tmp.vector)
            res <- rbind(res, c(species.names[species.index], final.tally, length(tmp.vector)))
        }
    }
    final.res <- data.frame(species=res[,1], region=res[,2], tot.records=res[,3])
    print(paste("temp", length(which(final.res$region==0))))
    print(paste("trop", length(which(final.res$region==1))))
    print(paste("both", length(which(final.res$region>0 & final.res$region<1))))
    return(final.res)
}

ecoregions.set <- fread("campBIOME30", sep=",")
latlongs <- fread("campGBIFRecordsCleanedCHECK")
tots <- cbind(ecoregions.set, latlongs$V3)
pp <- SummarizeBiomes(tots)
write.table(pp, file="camp.trop.temp", quote=FALSE, sep="\t", row.names=FALSE)
save(pp, file="biome.rough.Rsave")

####BIOME KEY####
#1 = Tropical and subtropical moist broadleaf forests
#2 = Tropical and subtropical dry broadleaf forests
#3 = Tropical and subtropical coniferous forests
#4 = Temperate broadleaf and mixed forests
#5 = Temperate coniferous forests
#6 = Boreal forests/Taiga
#7 = Tropical and subtropical grasslands, savannahs, and shrublands
#8 = Temperate grasslands, savannahs, and shrublands
#9 = Flooded grasslands and savannahs
#10 = Montane grasslands and shrublands
#11 = Tundra
#12 = Mediterranean forests, woodlands, and shrublands
#13 = Deserts and xeric shrublands
#14 = Mangroves
#98 = Lakes
#99 = Rock and ice

