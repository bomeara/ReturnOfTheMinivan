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

#' Download all no looping, new path to test
#'
#' They are all saved into data/pdfcache
#'
#' @param references.df data.frame from GetAllReferences
#' @return list of information
CacheAllPDFsImmediatelyDebug <- function(references.df) {
  fulltext::cache_options_set(full_path="/Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcachedebug")

  cache_all <- fulltext::ft_get(rev(references.df$DI))
  return(cache_all)
}


#' Debugging code to help Scott
#' @param cache_all output of CacheAllPDFsImmediately
#' @return DOIs that are failing
FindMissingDOIs <- function(cache_all) {
  paths <- cache_all$wiley$data$path
  full_path="/Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcache"
  baduns <- c()
  for (i in sequence(length(paths))) {
    if(is.null(paths[[i]]$type)) {
      output_file <- paste0(full_path, "/", gsub("\\.", "_", gsub('/', "_", paths[[i]]$id)), ".pdf")
      if(!file.exists(output_file)) {
        baduns <- append(baduns, paths[[i]]$id)

        fulltext::cache_options_set(full_path="/Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcachedebug")
        print(paths[[i]]$id)
        cache_all_bad <- fulltext::ft_get(paths[[i]]$id)
      }
    }
  }
  return(baduns)
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

#' Aggregate taxonomy by rank
#'
#' @param classifications Output of GetAllTaxonomy
#' @return data.frame with colnames different ranks
AggregateTaxonomy <- function(all_taxonomy) {
  results.df <- data.frame()
  for (taxon.index in sequence(nrow(all_taxonomy$resolved))) {
    ranks <- strsplit(all_taxonomy$resolved$classification_path_ranks[taxon.index], split="\\|")[[1]]
    elements <- strsplit(all_taxonomy$resolved$classification_path[taxon.index], split="\\|")[[1]]
    names(elements) <- ranks
    local.df <- as.data.frame(t(elements), stringsAsFactors=FALSE)
    results.df <- plyr::rbind.fill(results.df, local.df)
    if(taxon.index %% 100 == 0) {
      print(paste("done", taxon.index, "of", nrow(all_taxonomy$resolved)))
    }
  }
  return(results.df)
}

#' Get all plants from Catalogue of Life
#' @return data.frame of all plants
GetAllPlants <- function(extant_only=TRUE) {
  # Plantae is 208cf441fe2e1662376a9ce9e80782e1
  families <- taxize::downstream("208cf441fe2e1662376a9ce9e80782e1", db="col", downto="family", intermediate=FALSE)[[1]]
  families <- families[-which(grepl("Not assigned", families$childtaxa_name)),]
  all_species <- data.frame()
  for (family.index in sequence(nrow(families))) {
    print(paste("Doing family", family.index, families$childtaxa_name[family.index], "of", nrow(families)))
    genera <- taxize::downstream(families$childtaxa_id[family.index], db="col", downto="genus", intermediate=FALSE)[[1]]
    for (genus.index in sequence(nrow(genera))) {
      species <- taxize::downstream(genera$childtaxa_id[genus.index], db="col", downto="species", intermediate=FALSE)[[1]]
      if(any(species$childtaxa_extinct) & extant_only) {
        species <- species[!species$childtaxa_extinct,]
      }
      if(nrow(species) > 0) {
        all_species <- rbind(all_species, data.frame(family=families$childtaxa_name[family.index], genus=genera$childtaxa_name[genus.index], species=species$childtaxa_name, stringsAsFactors=FALSE))
      }
    }
  }
  return(all_species)
}

#' Count number of species in each genus, collapse GetAllPlants result to genus level
#'
#' Ranks are awful. But awfully useful, sadly.
#'
#' @param all_plants Output of GetAllPlants()
#' @return dataframe with columns for family, genus, and number of species
CountSpeciesByGenus <- function(all_plants) {
  final.df <- data.frame()
  unique_genera <- unique(all_plants$genus)
  for (genus.index in seq_along(unique_genera)) { # we're assuming no intra-plant genus homonyms. Code, don't fail us now! [actually, we check below]
    local.df <- subset(all_plants, all_plants$genus==unique_genera[genus.index])
    if(length(unique(local.df$family))>1) {
      print("Problem with this genus: in two families")
      print(local.df)
    } else {
      final.df <- rbind(final.df, data.frame(family=local.df$family[1], genus=local.df$genus[1], nspecies=nrow(local.df), stringsAsFactors=FALSE))
    }
  }
  return(final.df)
}

#' Get presence of matches
#'
#' Takes the described genera, sees if they're present or absent in another set of genera
#'
#' @param species_by_genus Output of CountSpeciesByGenus
#' @param genera genera to match, such as taxonomy_aggregated$genera
#' @param colname name for column of matches
GetMatches <- function(species_by_genus, genera, colname) {
  species_by_genus[,colname] <- species_by_genus$genus %in% genera
  return(species_by_genus)
}


#' Smith & Brown
#'
#' Uses the V0.1 release https://github.com/FePhyFoFum/big_seed_plant_trees/releases
#' From https://bsapubs.onlinelibrary.wiley.com/doi/full/10.1002/ajb2.1019
#'
#' We pull in the tree from genbank data only (GBMB.tre), resolve names.
ProcessSmithBrown <- function() {
  phy <- ape::read.tree("data/v0.1/GBMB.tre")
  return(AggregateTaxonomy(GetAllTaxonomy(data.frame(taxa=phy$tip.label, stringsAsFactors=FALSE))))
}

#' CalculateCoverage
#'
#' What fraction of genera or families are covered, and what proportion of total species each of these represents
#' @param input Result of GetMatches()
#' @return vector with the coverage
CalculateCoverage <- function(input) {
  result <- c()
  other_categories <- colnames(input)[-c(1:3)]
  total_species <- sum(input$nspecies)
  total_genera <- length(unique(input$genus))
  total_families <- length(unique(input$family))
  for (category.index in seq_along(other_categories)) {
    local <- input[input[,other_categories[category.index]],]
    local_total_species <- sum(local$nspecies)
    local_total_genera <- length(unique(local$genus))
    local_total_families <- length(unique(local$family))
    result <- c(result, paste(other_categories[category.index], "had genera with", local_total_species, "of", total_species, "so ", round(100*local_total_species/total_species,2), "percent, with", local_total_genera, "of", total_genera, "genera and ", local_total_families, "of", total_families, "families"))
  }
  print(result)
  return(result)
}

#' Get the genus only
#' @param A species name with an underscore between genus and species
GetGenus <- function(x) {
  return(strsplit(x, "_")[[1]][1])
}

#' Make tree of genera
#'
#' Using Smith and Brown's ALLMB.tre, merge species to genera
#' @return A chronogram of one species per genus
GenusTree <- function() {
  phy <- ape::read.tree("data/v0.1/GBMB.tre")
  taxa <- phy$tip.label
  genera <- unique(unname(sapply(taxa, GetGenus)))
  for (genus.index in seq_along(genera)) {
    to_kill <- taxa[grepl(paste0(genera[genus.index],"_"), taxa)][-1] #delete all but one
    if(length(to_kill)>0) {
      phy <- ape::drop.tip(phy, tip=to_kill)
    }
  }
  return(phy)
}

#' TNRS genus tree
#'
#' @param phy Chronogram of species, one per genus
#' @return Chronogram of genera
TNRSGenusTree <- function(phy) {
  all.tips <- unique(unname(sapply(phy$tip.label, GetGenus)))
  resolved_taxa <- data.frame()
  chunk_size <- 100
  for (i in seq(1, length(all.tips), chunk_size)) {

    seq_size <- chunk_size
    if ((i + seq_size) > length(all.tips)) seq_size <- length(all.tips) - i + 1

    resolved_taxa <- rbind(resolved_taxa, taxize::gnr_resolve(names=all.tips[i:i+seq_size], source=1, with_context=TRUE, best_match_only=TRUE)) #CoL
  }
  for (i in seq_along(phy$tip.label)) {
    if(phy$tip.label[i] %in% resolved_taxa$user_supplied_name) {
      phy$tip.label[i] <- resolved_taxa$matched_name[which(resolved_taxa$user_supplied_name==phy$tip.label[i])[1]]
    }
  }
  return(phy)
}

#' TNRS only genera
#' @param phy chronogram with species name already tnrs'ed
#' @return Chronogram of actual genera
TNRSGenusOnlyTree <- function(phy) {
  for (i in seq_along(phy$tip.label)) {
    phy$tip.label[i] <- GetGenus(phy$tip.label[i])
  }
  return(phy)
}

#' Get tip states
#'
#' Gets presence absence of a group in the AJB.
#' genus_binary is 0 or 1 (1=present)
#' family_three is 0=absent, 1=family present, not this genus, 2=family and genus present
#' genus_nspecies is number of species in the genus
#' family_nspecies is number of species in the family containing this genus
MakeTipData <- function(phy, species_by_genus_ajb_sb) {
  tip.data <- data.frame(taxon=phy$tip.label, family=NA, genus_binary=NA, family_three=NA, genus_nspecies=NA, family_nspecies=NA, stringsAsFactors=FALSE)
  for (i in seq_along(phy$tip.label)) {
    relevant.row <- which(species_by_genus_ajb_sb$genus==phy$tip.label[i])
    if(length(relevant.row)==1) {
      tip.data$genus_binary[i] <- ifelse(species_by_genus_ajb_sb$ajb[relevant.row],1,0)
      tip.data$genus_nspecies[i] <- species_by_genus_ajb_sb$nspecies[relevant.row]
      tip.data$family[i] <- species_by_genus_ajb_sb$family[relevant.row]
      family_rows <- which(species_by_genus_ajb_sb$family==species_by_genus_ajb_sb$family[relevant.row])
      if(length(family_rows)>0) {
        tip.data$family_three[i] <- tip.data$genus_binary[i] + ifelse(any(species_by_genus_ajb_sb$ajb[family_rows]),1,0)
        tip.data$family_nspecies <- sum(species_by_genus_ajb_sb$nspecies[family_rows])
      }
    }
  }
  return(tip.data)
}

# get rid of text cols
MakeTipDataForAde <- function(tip_data) {
  rownames(tip_data) <- tip_data$taxon
  tip_data <- tip_data[,-c(1,2)]
  return(tip_data)
}

PruneAll <- function(phy, tip_data_ade) {
  pruned <- geiger::treedata(phy=phy, data=tip_data_ade, sort=TRUE, warnings=FALSE)
  #pruned$phy <- ape::multi2di(pruned$phy)
  return(pruned)
}

PlotBullseye <- function(pruned, outfile="bull.pdf") {
  pdf(file=outfile, width=20, height=20)
  adephylo::bullseye(phy=pruned$phy, traits=pruned$data, legend=FALSE, alpha=0.5, cex=0.3)
  dev.off()
}

AceRecon <- function(pruned) {
  rates <- matrix(c(0, 1, 0, 0, 0, 1, 0, 0, 0), nrow=3, byrow=TRUE)
  tips <- as.numeric(pruned$data[,"family_three"])
  recon <- ape::ace(tips, multi2di(pruned$phy), type="discrete", model=rates)
  return(recon)
}
