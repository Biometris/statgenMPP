## New release

- Initial CRAN release

----

## Test environments

* local Windows 10 install, R 4.1.3
* winbuilder (develop)
* macbuilder
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results
There were no ERRORs or WARNINGs.

There was 1 NOTE

* checking CRAN incoming feasibility ... NOTE
Maintainer: 'Bart-Jan van Rossum <bart-jan.vanrossum@wur.nl>'

New submission

Possibly misspelled words in DESCRIPTION:
  IBD (23:70)
  MPP (23:44)
  QTL (2:8, 25:49)
  al (26:11)
  et (26:8)
  
These are spelled correctly.  

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1073/pnas.1100465108
    From: inst/doc/QTLMapping_in_MultiParentalPopulations.html
    Status: 503
    Message: Service Unavailable
  URL: https://www.jstor.org/stable/29713
    From: inst/doc/QTLMapping_in_MultiParentalPopulations.html
    Status: 403
    Message: Forbidden

Both urls work correctly when opened from my browser.
