#' @importFrom ggplot2 stat
NULL
#' Returns histogram for grouped data
#'
#' @param x Numeric vector of values that will be plotted.
#' @param group A vector that contains the grouping variable. It can be a 
#' numeric, factor, or character vector.
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
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' sce <- labelQC(sce)
#' cytofHist(scores(sce, 'bead'), label(sce))
#' 
#' @export
cytofHist <- function(x, group, type = c("count", "density"), na.rm = FALSE, 
                      title = NULL) {
    
    if (!methods::is(x, "numeric")) {
        stop("x must be a numeric vector")
    }
    
    type <- match.arg(tolower(type), choices = c("count", "density"))
    
    tmp <- data.frame(x, group = group)
    cols <- c("#69b3a2", "#404080", "#80b1d3", "#d4b9da", "#fdd0a2")
    
    if (na.rm) {
        tmp <- tmp[!is.na(tmp$group), ]
    }
    
    cols <- cols[seq_along(unique(tmp$group))]
    
    g <- ggplot2::ggplot(tmp, ggplot2::aes(x = x, fill = factor(group)))
    
    if (type == "count") {
        g <- g + ggplot2::geom_histogram(color = "#e9ecef", alpha=0.6,
                                         position = 'identity', bins = 100)
    } else if (type == "density") {
        density <- NULL
        g <- g + ggplot2::geom_histogram(color="#e9ecef", alpha=0.6,
                                         position = 'identity',
                                         bins = 100, 
                                         ggplot2::aes(y = stat(density)))
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
