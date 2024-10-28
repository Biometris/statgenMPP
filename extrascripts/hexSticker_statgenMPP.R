# note: run usethis::use_logo() to generate the generate correct size

library(hexSticker)
library(ggplot2)
library(LMMsolver)
library(spam)
library(fields)
library(maps)
library(dplyr)
library(gridExtra)

markerFiles <- c(system.file("extdata/multipop", "AxB.txt",
                             package = "statgenMPP"),
                 system.file("extdata/multipop", "AxC.txt",
                             package = "statgenMPP"))

## Specify file containing map.
# Both crosses use the same map file.
mapFile <- system.file("extdata/multipop", "mapfile.txt",
                       package = "statgenMPP")

## Read phenotypic data
phenoDat <- read.delim(system.file("extdata/multipop", "AxBxCpheno.txt",
                                   package = "statgenMPP"))
# Check contents.
head(phenoDat)

## Perform IBD calculations.
ABCMPP <- calcIBDMPP(crossNames = c("AxB", "AxC"),
                     markerFiles = markerFiles,
                     pheno = phenoDat,
                     popType = "F4DH",
                     mapFile = mapFile,
                     evalDist = 5)

ABCMQM <- selQTLMPP(MPPobj = ABCMPP,
                    trait = "yield",
                    threshold = 3)

#p <- plot(ABCMQM, plotType = "QTLProfileExt")
p <- plot(ABCMQM, plotType = "QTLProfile", chr=c(1,2,3))
#p <- plot(maizeSQM, plotType = "QTLProfile",chr=c(1,2,3,4))
#p <- p + theme_void() + theme_transparent()
#p <- plot(ABCMQM, plotType = "parEffs")
p <- p + geom_text(x=80, y=18, label="QTL 1", size=8, col='darkblue')
p <- p + geom_text(x=180, y=11.2, label="QTL 2", size=8, col='darkblue')
p <- p + geom_text(x=260, y=14.2, label="QTL 3", size=8, col='darkblue')

p <- p + theme_transparent() + theme(axis.title=element_blank())
p

#sticker(p, package="statgenMPP", p_size=20, s_x=1.0, s_y=.85, s_width=1.2, s_height=0.9,
#        h_fill = "lightyellow", p_color="blue",h_color="green",
#        filename="statgenMPP_hexSticker.png")

sticker(p, package="statgenMPP", p_size=20, s_x=0.87, s_y=.7, s_width=1.4, s_height=1.0,
        h_fill = "lightyellow", p_color="blue",h_color="green",
        filename="statgenMPP_hexSticker.png")

# sticker(p, package="statgenMPP", p_size=18, s_x=0.9, s_y=.8, s_width=1.4, s_height=1.0,
#         h_fill = "lightyellow", p_color="blue",h_color="green",
#         filename="statgenMPP_hexSticker.png")




