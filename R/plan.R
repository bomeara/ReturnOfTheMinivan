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
  save_all_taxonomy_csv = write.csv(all_taxonomy$resolved, file=file_out('all_taxonomy.csv'))
)
