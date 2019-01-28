
########FIGURE MAKER#########

#written by Jeremy M. Beaulieu

# library(ape)
# library(corHMM)
# library(RColorBrewer)
# source("plot-alt.R")
# source("plot-alt-util.R")
# source("plot-alt-extra.R")


MakePlot<-function(phy, data, pars, col.vec){

    tmp <- rayDISC(phy, data, p=pars, model="ARD", node.states="marginal")

    phy <- tmp$phy
    phy<- ladderize(phy, right=FALSE)
    data <- data.frame(data[,2], data[,2], row.names=data[,1])
	data <- data[phy$tip.label,]
	nb.tip <- length(phy$tip.label)
	nb.node <- phy$Nnode

    tips <- as.numeric(data[,1])+1
    comp <- numeric(Nedge(phy))
    comp[match(1:(Ntip(phy)), phy$edge[,2])] <- tips
    comp[match((2+Ntip(phy)):(Nedge(phy)+1), phy$edge[,2])] <- as.numeric(phy$node.label[-1])
	colors<-col.vec

	#pdf(file=file.name, width=8, height=10)
	plot2.phylo(phy,show.tip.label=TRUE,lend=1, edge.width=1, edge.color=colors[comp], type="fan", cex=0.05)
#	dev.off()
}



GetColor <- function(x) {
	if(which.max(x)==1) {
		return("gray")
	}
	normalizing <- x[2] + x[4]
	weight.blue <- x[2]/normalizing
	weight.red <- x[4]/normalizing
	return(rgb(weight.red, 0, weight.blue))
}


MakePlotUncertain <- function(hisse.states.obj, file.name="plotUncertainHiSSE.pdf" ){
	phy<-ladderize(hisse.states.obj$phy, right=FALSE)
	node.mat <- hisse.states.obj$node.mat
	tip.mat <- hisse.states.obj$tip.mat
	tot.mat <- rbind(tip.mat, node.mat[-1,])
	comp <- matrix(0, dim(tot.mat)[1], dim(node.mat)[2])
	comp[match(1:(Ntip(phy)), phy$edge[,2]),] <- tip.mat
	comp[match((2+Ntip(phy)):(Nedge(phy)+1), phy$edge[,2]),] <- node.mat[-1,]
	pdf(file=file.name, width=10, height=10)
	plot2.phylo(phy,show.tip.label=FALSE,lend=1, edge.color=apply(comp, 1, GetColor), type="fan")
	dev.off()
}


### Example ###
##################################################################
#data(primates)
#pars <- c(0.01955962, 0.005879677)
#MakePlot(primates$tree, primates$trait, pars=pars, #col.vec=c("red", "blue"))
