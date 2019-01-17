my_plan <- drake_plan(
  references.df = GetAllReferences(),
  references.txt = DownloadAndExtractAllPDFs(references.df)
)

my_plan_immediately <- drake_plan(
  references.df = GetAllReferences(),
  cache_all = CacheAllPDFsImmediately(references.df),
  all_pdfs = CacheRemainingPDFs(cache_all),
  all_names = GetAllNames(all_pdfs)
)
