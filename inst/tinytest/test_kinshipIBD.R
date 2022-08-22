### Test kinshipIBD

set.seed(1234)
## Create marker data.
markers <- array(runif(45), dim = c(3, 5, 3),
                 dimnames = list(paste0("G", 1:3), paste0("M", 1:5),
                                 paste0("P", 1:3)))
markers <- simplify2array(apply(X = markers, MARGIN = 2, FUN = function(x) {
  x / rowSums(x)
}, simplify = FALSE))
markers <- aperm(markers, c(1, 3, 2))

## Construct map.
map <- data.frame(chr = c(1, 1, 2, 2, 2), pos = 1:5,
                  row.names = paste0("M", 1:5))


## General checks

expect_error(kinshipIBD(markers = 1),
             "markers should be a 3 dimensional array")
expect_error(kinshipIBD(markers = markers, map = 1),
             "map should be a data.frame")
expect_error(kinshipIBD(markers = markers, map = map[, 1, drop = FALSE]),
             "chr and pos should be columns in map")

## Chromosome specific
KChrSpec <- kinshipIBD(markers = markers, map = map)

expect_inherits(KChrSpec, "list")
expect_equal(length(KChrSpec), 2)
expect_equal(names(KChrSpec), c("1", "2"))
expect_inherits(KChrSpec[[1]], "matrix")
expect_equal(dim(KChrSpec[[1]]), c(3, 3))

expect_equal(as.numeric(KChrSpec[[1]]),
             c(0.450980778777623, 0.318347657907409, 0.274145151365274,
               0.318347657907409, 0.342980515400794, 0.343517588921626,
               0.274145151365274, 0.343517588921626, 0.423994558163868))

## Non chromosome specific

K <- kinshipIBD(markers = markers, chrSpecific = FALSE)

expect_inherits(K, "matrix")
expect_equal(dim(K), c(3, 3))

expect_equal(as.numeric(K),
             c(0.442181446699766, 0.323298364111491, 0.291059186451823,
               0.323298364111491, 0.386092061743138, 0.353213553269824,
               0.291059186451823, 0.353213553269824, 0.404976882045675))
