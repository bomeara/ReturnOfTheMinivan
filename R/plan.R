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
  coverage = CalculateCoverage(species_by_genus_ajb_sb),
  genus_tree = GenusTree(),
  tnrs_genus_tree = TNRSGenusTree(genus_tree),
  save_tnrs_genus_tree = ape::write.tree(tnrs_genus_tree, file=file_out("tnrs_genus_tree.tre")),
  tnrs_genus_only_tree = TNRSGenusOnlyTree(tnrs_genus_tree),
  save_tnrs_genus_only_tree = ape::write.tree(tnrs_genus_only_tree, file=file_out("tnrs_genus_only_tree.tre")),
  tip_data=MakeTipData(tnrs_genus_only_tree, species_by_genus_ajb_sb),
  tip_data_ade = MakeTipDataForAde(tip_data),
  pruned = PruneAll(tnrs_genus_only_tree, tip_data_ade),
  bullseye_plot = PlotBullseye(pruned, outfile=file_out("bull.pdf")),
  genus_recon = AceRecon2(pruned),
  branch_plot = PlotRecon(pruned,genus_recon, outfile=file_out("recon.pdf")),
  jeremy_plot = PlotJeremy(pruned, outfile=file_out("jeremy.pdf"))
)
