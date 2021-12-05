#' Returns histogram for grouped data
#'
#' @param x The vector of values that will be plotted.
#' @param group A vector that contains the grouping variable.
#' @param type Either "count" or "density". The "count" selection
#' keeps the groups on the same scale. The "density" option will
#' over emphasize the group with the fewest observations. This
#' is helpful when identifying where certain subgroups are relative
#' to the majority of the data.
#' @param na.rm TRUE if NAs should be removed prior to plotting.
#' FALSE if they should remain. The NAs will be plotted as a
#' separate group if they are not removed.
#' @param title Optional title for the plot
#'
#' @return A \code{ggplot2} histogram.
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' nn <- NN(x)
#' labels <- qcDataFrame(x)
#' beads <- initialBead(x, labels = labels)
#' cytofHist(beads$Bead1, beads$init)
#'
#' @export
cytofHist <- function(x, group, type = "count", na.rm = FALSE, title = NULL) {

  tmp <- data.frame(x, group = group)
  cols <- c("#69b3a2", "#404080", "#80b1d3", "#d4b9da", "#fdd0a2")

  if (na.rm) {
    tmp <- tmp[!is.na(tmp$group), ]
  }

  cols <- cols[1:length(unique(tmp$group))]

  g <- ggplot2::ggplot(tmp, ggplot2::aes(x = x, fill = factor(group)))

  if (type == "count") {
    g <- g + ggplot2::geom_histogram(color = "#e9ecef", alpha=0.6,
                                     position = 'identity', bins = 100)
  } else if (type == "density") {
    g <- g + ggplot2::geom_histogram(color="#e9ecef", alpha=0.6,
                                     position = 'identity',
                                     bins = 100, ggplot2::aes(y = ..density..))
  } else {
    warning("invalid type - using counts for histogram")
    g <- g + ggplot2::geom_histogram(color = "#e9ecef", alpha=0.6,
                                     position = 'identity', bins = 100)
  }

  g <- g + ggplot2::scale_fill_manual(values = cols) +
    hrbrthemes::theme_ipsum() +
    ggplot2::labs(fill = "") +
    ggplot2::theme_bw()

  if(!is.null(title)) {
    g <- g + ggplot2::ggtitle(title)
  }

  g

}
