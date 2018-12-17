# ReturnOfTheMinivan
Data for response to Edwards &amp; Donoghue's response to Beaulieu and O'Meara. We have more cupholders, and they come filled with data! [and cheerios]

data/AJB_recs1to500.txt and data/AJB_recs501to1000.txt have the info on the most recent 1000 articles in AJB (downloaded on Dec. 17 from Web of Knowledge. part of 2013 to Nov 2018.

* 877 are articles
* 90 are editorials
* 14 are reviews
* 9 are corrections
* 6 are letters
* 4 are news

Top authors: Soltis DE (15), Soltis PS (15), Ashman TL (10), Pires JC (9), Rothwell GW (9), Smith SA (8), Escudero M (8), Rudall PJ (7), Temescu AMF (7)


DI field has DOI

library(fulltext)
ft_get("10.1002/ajb2.1180") # gets the PDF
ft_extract(ft_get("10.1002/ajb2.1180"))$wiley$data # is all the text

You can edit this to get rid of the references field, then pass to rphylotastic to get names
