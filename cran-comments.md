## Second release

- Improvements in the algorithm to increase performance and reduce memory allocation. Also allow for parallel computing.

----

## Test environments

* local Windows 10 install, R 4.2.1
* winbuilder (develop)
* macbuilder
* Ubuntu (on github actions, devel and release)
* macOS (on github actions, release)
* R-hub (devel and release)

----

## R CMD check results

There were no ERRORs or WARNINGs.

There was 1 NOTE

Found the following (possibly) invalid URLs:
  URL: https://doi.org/10.1073/pnas.1100465108
    From: inst/doc/QTLMapping_in_MultiParentPopulations.html
    Status: 503
    Message: Service Unavailable
  URL: https://www.jstor.org/stable/29713
    From: inst/doc/QTLMapping_in_MultiParentPopulations.html
    Status: 403
    Message: Forbidden
  
These links work correctly when opened from my browser.
