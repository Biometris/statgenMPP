## Check with winbuilder - develop only.
devtools::check_win_devel()

## Check with macbuilder
devtools::check_mac_release()


rhub::platforms()

## Check with rhub.
rhub::check_for_cran(path = "C:/Projects/R_packages/statgenMPP/")

## Check on rhub for roldrel. This crashed for v.1.0.0
rhub::check_with_roldrel(path = "C:/Projects/R_packages/statgenMPP/")

## Rebuild readme.
devtools::build_readme()

## Build site for local check.
pkgdown::clean_site()
pkgdown::build_site()

## Submit to CRAN
devtools::release()

## Code coverage - local.
detach("package:statgenMPP", unload = TRUE)
covr::gitlab()
library(statgenMPP)


