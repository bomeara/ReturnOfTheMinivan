devtools::install_github("thej022214/corHMM") #dev version of corHMM
library(ape)
library(corHMM)
library(geiger)

Discretize <- function(x, threshold=0.5) {
  x[x<threshold] <- 0
  x[x>=threshold] <- 1
  return(x)
}

DoSingleAnalysis <- function(tree, data, threshold=0.5, ...) {
   data.local <- data
   data.local$region <- Discretize(data.local$region, threshold)
   data.local <- data.local[,-3]
   result <- corHMM(tree, data.local, root.p="maddfitz", ...)
   return(result)
}
