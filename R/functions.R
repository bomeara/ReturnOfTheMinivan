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
  fulltext::cache_options_set(full_path="/Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcache")

  cache_all <- fulltext::ft_get(rev(references.df$DI))
  return(cache_all)
}

#' Hack to get remaining ones
#'
#' Wiley sometimes wants to give ePDFs. Curse them.
#'
#' Thanks to Scott Chamberlain for the crul workaround
#'
#' @param cache_all output of CacheAllPDFsImmediately
#' @return Paths to all PDFs
CacheRemainingPDFs <- function(cache_all) {
  paths <- cache_all$wiley$data$path
  full_path="/Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcache"
  for (i in sequence(length(paths))) {
    if(is.null(paths[[i]]$type)) {
      output_file <- paste0(full_path, "/", gsub("\\.", "_", gsub('/', "_", paths[[i]]$id)), ".pdf")
      if(!file.exists(output_file)) {
        print(output_file)
        cli = crul::HttpClient$new(url=paste0("https://bsapubs.onlinelibrary.wiley.com/doi/pdf/", paths[[i]]$id), headers = list(Accept = "application/pdf", "CR-Clickthrough-Client-Token" =  Sys.getenv("CROSSREF_TDM")), opts = list(followlocation=1))
        try(res <- cli$get(disk=output_file))
      }
    }
  }
  return(list.files(path=full_path, full.names=TRUE))
}


#' Extract scientific names from files
#'
#' @param file_paths All paths to the PDFs
#' @return list of taxon names (taxa) and the files they came from (files)
GetAllNames <- function(file_paths) {
  #file_paths <- file_paths[1:3]
  all_names <- vector("list", length(file_paths))
  all_files <- vector("list", length(file_paths))
  good_ones <- rep(FALSE, length(file_paths))
  for (i in seq_along(file_paths)) {
    taxon_names <- NULL
    try(taxon_names <- rphylotastic::file_get_scientific_names(file_paths[i]))
    if(!is.null(taxon_names)) {
      all_names[[i]] <- taxon_names
      all_files[[i]] <- file_paths[i]
      good_ones[i] <- TRUE
    }
  }
  return(list(taxa=all_names[good_ones], files=all_files[good_ones]))
}

#' Get the taxonomy for all the taxa in papers
#'
#' @param all_names Result from GetAllNames
#' @return output of datelife::classification_paths_from_taxonomy()
GetAllTaxonomy <- function(all_names) {
  unique_taxa <- unique(unlist(all_names$taxa))
  classifications <- datelife::classification_paths_from_taxonomy(unique_taxa)
  return(classifications)
}
