my_plan <- drake_plan(
  references.df = GetAllReferences(),
  references.txt = DownloadAndExtractAllPDFs(references.df)
)
