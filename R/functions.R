#' Load all the references, format into a data.frame
#'
#' @return a dataframe with all the references
GetAllReferences <- function() {
  setwd("data")
  files <- list.files(pattern=".txt")
  references.df <- data.frame()
  for (i in seq_along(files)) {
    local.files <- bibliometrix::readFiles(files[i])
    local.files <- bibliometrix::convert2df(local.files, dbsource = "isi", format = "plaintext")
    if(i==1) {
      references.df <- local.files
    } else {
      references.df <- plyr::rbind.fill(references.df, local.files)
    }
  }
  setwd("..")
  return(references.df)
}

#' Download all PDFs
#'
#' They are all saved into data/pdfcache
#'
#' @param references.df data.frame from GetAllReferences
#' @return list of information
DownloadAndExtractAllPDFs <- function(references.df) {
  fulltext::cache_options_set(path="../../../../../Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcache") #yes, this is stupid
  references.txt <- vector("list", nrow(references.df))
  for (i in seq_along(references.df$DI)) {
    try(references.txt[[i]] <- fulltext::ft_extract(ft_get(references.df[i,"DI"])))
    Sys.sleep(30+runif(1,0,80))
  }
  return(references.txt)
}

#' Download all no looping
#'
#' They are all saved into data/pdfcache
#'
#' @param references.df data.frame from GetAllReferences
#' @return list of information
CacheAllPDFsImmediately <- function(references.df) {
  fulltext::cache_options_set(path="../../../../../Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcache") #yes, this is stupid
  cache_all <- fulltext::ft_get(references.df$DI[1:3])
  return(cache_all)
}
