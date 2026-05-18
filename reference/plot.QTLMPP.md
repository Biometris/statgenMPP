# Plot function for the class `QTLMPP`

Creates a plot of an object of S3 class `QTLMPP`. The following types of
plot can be made:

- QTLProfile:

  A QTL profile plot, i.e. a plot of \\-log10(p)\\ values per marker.

- parEffs:

  A plot of effect sizes and directions per parent.

- QTLRegion:

  A plot highlighting the QTLs found on the genetic map.

- QTLProfileExt:

  A combination of the QTL profile and parental effects plotted above
  each other.

- parCIs:

  A plot of the effect sizes and confidence intervals for each parent
  for each of the QTLs found.

See details for a detailed description of the plots and the plot options
specific to the different plots.

## Usage

``` r
# S3 method for class 'QTLMPP'
plot(
  x,
  ...,
  plotType = c("QTLProfile", "parEffs", "QTLRegion", "QTLProfileExt", "parCIs"),
  title = NULL,
  output = TRUE
)
```

## Arguments

- x:

  An object of class `QTLMPP`.

- ...:

  further arguments to be passed on to the actual plotting functions.

- plotType:

  A character string indicating the type of plot to be made. One of
  "QTLProfile", "parEffs", and "QTLRegion".

- title:

  A character string, the title of the plot.

- output:

  Should the plot be output to the current device? If `FALSE`, only a
  ggplot object is invisibly returned.

## QTL Profile plot

A profile of all marker positions and corresponding \\-log10(p)\\ values
is plotted. QTLs found are highlighted with red dots. The threshold is
plotted as a horizontal line. If there are previously known marker
effects, false positives and true negatives can also be marked.  
Extra parameter options:

- `xLab`:

  A character string, the x-axis label. Default = `"Chromosomes"`

- `yLab`:

  A character string, the y-axis label. Default = `-log10(p)`

- `effects`:

  A character vector, indicating which markers correspond to a real
  (known) effect. Used for determining true/false positives and false
  negatives. True positives are colored green, false positives orange
  and false negatives yellow.

- `colPalette`:

  A color palette used for plotting. Default coloring is done by
  chromosome, using black and grey.

- `signLwd`:

  A numerical value giving the thickness of the points that are
  false/true positives/negatives. Default = 0.6

- `chr`:

  A vector of chromosomes to be plotted. By default, all chromosomes are
  plotted. Using this option allows restricting the plot to a subset of
  chromosomes.

## Parental effects Plot

A plot of effect sizes for each of the parents for the QTLs found is
created.  
Extra parameter options:

- `xLab`:

  A character string, the x-axis label. Default = `"Chromosomes"`

- `yLab`:

  A character string, the y-axis label. Default = `"Parents"`

- `chr`:

  A vector of chromosomes to be plotted. By default, all chromosomes are
  plotted. Using this option allows restricting the plot to a subset of
  chromosomes.

## QTL Region Plot

A plot of effect is created on which the positions of the QTLs found in
the QTL detection are highlighted in red.  
No extra parameter options.

## Extended QTL Profile Plot

An extended version of the QTL Profile Plot, in which the QLT profile
plot is combined with the parental effect plot to make it easier to
assess the effects for each specific QTL found.  
Extra parameter options:

- `chr`:

  A vector of chromosomes to be plotted. By default, all chromosomes are
  plotted. Using this option allows restricting the plot to a subset of
  chromosomes.

## Parental Confidence Interval Plot

A plot of the parental effects and the confidence intervals around them
for each of the QTLs found. Positive significant effects are shown in
red, negative significant effects in blue, and all other effects are
shown in black.

## Examples

``` r
if (FALSE) { # \dontrun{
## Read phenotypic data.
pheno <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                               package = "statgenMPP"))
## Rename first column to genotype.
colnames(pheno)[1] <- "genotype"


## Compute IBD probabilities for simulated population - AxB, AxC.
ABC <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                  markerFiles = c(system.file("extdata/multipop", "AxB.txt",
                                              package = "statgenMPP"),
                                  system.file("extdata/multipop", "AxC.txt",
                                              package = "statgenMPP")),
                  pheno = pheno,
                  popType = "F4DH",
                  mapFile = system.file("extdata/multipop", "mapfile.txt",
                                        package = "statgenMPP"),
                  evalDist = 5)

## Multi-QTL Mapping.
ABC_MQM <- selQTLMPP(ABC, trait = "yield")

## QTL Profile plot.
plot(ABC_MQM, plotType = "QTLProfile")

## Plot of parental effects for QTLs found.
plot(ABC_MQM, plotType = "parEffs")

## Plot of genetic map highlighting positions of QTLs found.
plot(ABC_MQM, plotType = "QTLRegion")

## Extended QTL Profile plot.
plot(ABC_MQM, plotType = "QTLProfileExt")

## Parental confidence interval plot.
plot(ABC_MQM, plotType = "parCIs")
} # }
```
