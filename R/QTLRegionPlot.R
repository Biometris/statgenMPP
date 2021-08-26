QTLRegionPlot <- function(signSnp,
                          map,
                          chrBoundaries,
                          title = NULL,
                          output = TRUE) {
  signSnp[["notation"]] <- paste0(signSnp[["chr"]], "@", signSnp[["pos"]])

  df <- data.frame(x = unique(map[["chr"]]), y = 0,
                   xend = unique(map[["chr"]]), yend = chrBoundaries[["x"]])

  if (is.null(title)) {
    title <- "Genetic map"
  }
  p <- ggplot2::ggplot(map, ggplot2::aes(x = chr, y = pos)) +
    ggplot2::geom_segment(data = df, ggplot2::aes(x = x, y = y,
                                                xend = xend, yend = yend),
                          lineend = "round", size = 8, color = "grey90") +
    ggplot2::geom_segment(data = df, ggplot2::aes(x = x, y = y,
                                                  xend = xend, yend = yend)) +
    ggplot2::geom_point(pch = "_", size = 4) +
    ggplot2::geom_point(data = signSnp, color = "red", pch = "_", size = 5) +
    ggplot2::geom_segment(data = signSnp,
                          ggplot2::aes(x = chr + 0.25, y = pos,
                                       xend = chr + 0.02, yend = pos),
                          arrow = ggplot2::arrow(length = ggplot2::unit(0.3, "cm")),
                          color = "red") +
    ggplot2::annotate(geom="text", x = signSnp$chr + 0.5,
                      y = signSnp$pos, label = signSnp$notation,
                      color = "red", size = 2) +
    ggplot2::labs(title = title, x = "Chromosome", y = "Position") +
    ggplot2::scale_x_discrete(limits = factor(unique(map$chr)),
                              labels = as.character(unique(map$chr))) +
    ggplot2::scale_y_reverse() +
    ggplot2::theme(panel.background = ggplot2::element_blank(),
                   plot.background = ggplot2::element_blank(),
                   axis.ticks.x = ggplot2::element_blank(),
                   strip.background = ggplot2::element_blank(),
                   panel.border = ggplot2::element_blank(),
                   axis.line.y = ggplot2::element_line(size = 0.5),
                   plot.title = ggplot2::element_text(hjust = 0.5))
  if (output) {
    plot(p)
  }
  invisible(p)
}
