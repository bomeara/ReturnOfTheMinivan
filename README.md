# ReturnOfTheMinivan
Data for response to Donoghue &amp; Edward's response to Beaulieu and O'Meara. We have more cupholders, and they come filled with data! [and cheerios]

data/AJB_recs1to500.txt and data/AJB_recs501to1000.txt have the info on the most recent 1000 articles in AJB (downloaded on Dec. 17 from Web of Knowledge. part of 2013 to Nov 2018.

* 877 are articles
* 90 are editorials
* 14 are reviews
* 9 are corrections
* 6 are letters
* 4 are news

Top authors: Soltis DE (15), Soltis PS (15), Ashman TL (10), Pires JC (9), Rothwell GW (9), Smith SA (8), Escudero M (8), Rudall PJ (7), Temescu AMF (7)


DI field has DOI

```
library(bibliometrix)
library(fulltext)


#setwd("/Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data")

files <- list.files(pattern=".txt")
all.files <- data.frame()
for (i in seq_along(files)) {
  local.files <- bibliometrix::readFiles(files[i])
  local.files <- bibliometrix::convert2df(local.files, dbsource = "isi", format = "plaintext")
  if(i==1) {
    all.files <- local.files
  } else {
    all.files <- plyr::rbind.fill(all.files, local.files)
  }
}


fulltext::cache_options_set(path="../../../../../Users/bomeara/Documents/MyDocuments/GitClones/ReturnOfTheMinivan/data/pdfcache")
#fulltext::ft_get("10.1002/ajb2.1180") # gets the PDF
text <- fulltext::ft_extract(ft_get(all.files[1,"DI"]))$wiley$data # is all the text
```

You can edit this to get rid of the references field, then pass to rphylotastic to get names

## To run

This uses [drake](https://ropenscilabs.github.io/drake-manual/index.html) to run the analyses. The necessary functions, packages, etc are described in scripts in the R directory. To rerun the analysis, just `source("make.R")` within R (or, better so it can run unattended: `nohup R CMD BATCH make.R > nohup.out &`).
