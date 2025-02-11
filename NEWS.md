# statgenMPP 1.0.4

* `readRABBITMPP` now reads the output files producted by open source version of RABBIT, available at <https://github.com/Biometris/RABBIT>. Older files created using Mathematica are still read as well.
* A new plot (`plotType = "parCIs"`) is now available for plotting the parental effects and the confidence intervals around them.

# statgenMPP 1.0.3

* The output of `selQTLMPP` now contains the standard errors of the effects. The standard errors are added in both the full output (GWAResult) and the table with QTLs (signSnp).
* The output of the two plots in QTLProfileExt is properly aligned again.
* Functions no longer rely on soft-deprecated ggplot2 functions.

# statgenMPP 1.0.2

* A minor bug in the summary for `gDataMPP` objects is fixed. Now the correct number of genotypes is shown also when there is only one cross.
* A citation file was added.

# statgenMPP 1.0.1.1

* Internal fix. No user visible changes.

# statgenMPP 1.0.1

* The logic in the scan process is improved, leading to much lower computation times. * Computations can now be done in parallel by registering a parallel back-end and setting `parallel = TRUE` in the `selQTLMPP` function.
* It is now possible to include a polygenic effect in the models to correct for background genetic information.
* A new function `createGDataMPP` is added to directly create a `gDataMPP` object from an `IBDprob` object (as created by `statgenIBD`) and a data.frame with phenotypic data.
* `readRABBIT` is renamed to `readRABBITMPP` to prevent future conflicts with a similar function in `statgenIBD`
* The trait column in the example data is renamed from "pheno" to "yield".

# statgenMPP 1.0.0

* Initial CRAN release
