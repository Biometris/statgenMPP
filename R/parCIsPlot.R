#' Create a plot with effect sizes
#'
#' Create a plot of QTL effect sizes and standard errors for each parent for
#' selected snps.
#'
#' @keywords internal
parCIsPlot <- function(parCIsDat,
                       xLab = "Effect",
                       yLab = "Parents",
                       title = NULL,
                       trait = NULL,
                       ...,
                       output = TRUE) {
  parCIsDat[["sign"]] <- as.factor(
    ifelse(parCIsDat[["effect"]] - 2 * parCIsDat[["seEffect"]] < 0 &
             parCIsDat[["effect"]] + 2 * parCIsDat[["seEffect"]] > 0, 0, 1))
  parCIsDat[["snp"]] <- factor(parCIsDat[["snp"]],
                               levels = unique(parCIsDat[["snp"]]))
  ## Create title.
  if (is.null(title)) {
    title <- paste("Parental effects at QTLs for", trait)
  }
  if (nrow(parCIsDat) > 0) {
    p <- ggplot2::ggplot(parCIsDat,
                         ggplot2::aes(x = .data[["effect"]],
                                      y = .data[["parent"]],
                                      color = .data[["sign"]])) +
      ggplot2::geom_point(color = "black") +
      ggplot2::geom_errorbarh(ggplot2::aes(xmax = .data[["effect"]] + 2 * .data[["seEffect"]],
                                           xmin = .data[["effect"]] - 2 * .data[["seEffect"]]),
                              height = 0.2) +
      ## Trick to assure symmetrical x-axis.
      ggplot2::geom_blank(ggplot2::aes(xmin = - .data[["effect"]] - 2 * .data[["seEffect"]],
                                       xmax = - .data[["effect"]] + 2 * .data[["seEffect"]])) +
      ggplot2::geom_vline(xintercept = 0, linetype = "dashed") +
      ggplot2::labs(x = xLab, y = yLab, title = title) +
      ggplot2::facet_wrap(ggplot2::vars(.data[["snp"]])) +
      ggplot2::scale_y_discrete(limits = rev) +
      ggplot2::scale_color_manual(values = c("0" = "black", "1" = "red"),
                                  guide = "none") +
      ggplot2::theme(panel.background = ggplot2::element_blank(),
                     plot.background = ggplot2::element_blank(),
                     strip.background = ggplot2::element_blank(),
                     panel.border = ggplot2::element_rect(fill = NA,
                                                          color = "black",
                                                          size = 0.5,
                                                          linetype = "solid"),
                     plot.title = ggplot2::element_text(hjust = 0.5))
  }
  if (output) {
    plot(p)
  }
  invisible(p)
}

