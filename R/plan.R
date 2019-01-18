pkgconfig::set_config("drake::strings_in_dots" = "literals")

my_plan <- drake_plan(
  references.df = GetAllReferences(),
  references.txt = DownloadAndExtractAllPDFs(references.df)
)

my_plan_immediately <- drake_plan(
  references.df = GetAllReferences(),
  cache_all = CacheAllPDFsImmediately(references.df),
  all_pdfs = CacheRemainingPDFs(cache_all),
  all_names = GetAllNames(all_pdfs),
  save_all_names = save(all_names, file=file_out('all_names.rda')),
  all_taxonomy = GetAllTaxonomy(all_names),
  save_all_taxonomy = save(all_taxonomy, file=file_out('all_taxonomy.rda')),
  save_all_taxonomy_csv = write.csv(all_taxonomy$resolved, file=file_out('all_taxonomy.csv')),
  taxonomy_aggregated = AggregateTaxonomy(all_taxonomy),
  save_aggregated = save(taxonomy_aggregated, file=file_out('taxonomy_aggregated.rda')),
  all_plants = GetAllPlants(),
  save_all_plants = save(all_plants, file=file_out("all_plants.rda")),
  species_by_genus = CountSpeciesByGenus(all_plants),  # note: five genera were in more than one family, and were tossed
  species_by_genus_ajb = GetMatches(species_by_genus, taxonomy_aggregated$genus, "ajb"),
  taxonomy_aggregated_smith_brown = ProcessSmithBrown(),
  species_by_genus_ajb_sb = GetMatches(species_by_genus_ajb, taxonomy_aggregated_smith_brown$genus, "smith_brown"),
  coverage = CalculateCoverage(species_by_genus_ajb_sb)
)
