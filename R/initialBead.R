#' Preliminary bead classification
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}}. 
#'
#' @return A \code{SingleCellExperiment} that contains the bead score and the 
#' bead designation for each event. This information is stored in the 
#' \code{score} and \code{initial} objects in the colData for the 
#' \code{SingleCellExperiment}. 
#'
#' @details
#' The beads are typically the first cell classification that is done. The
#' different event types are labeled iteratively so the \code{labels}
#' vector in the colData will contain all of the labels and 
#' probabilities computed up to this point. Only events that 
#' have a "cell" label can be assigned an initial event classification of
#' "bead". This function computes a score that assesses how much an event
#' looks like a bead and then fits a mixture model to assign each event 
#' a class of 1 for a bead, -1 for an event that is not a bead, or 0 
#' for undetermined or previously assigned to a different event type. 
#' The score is recorded in the \code{score} object in the colData and 
#' the initial classification is recorded in the \code{initial} part of 
#' the colData. 
#' 
#' Each bead channel should classify into two fairly clear groups where one
#' is the beads and the other is non-beads. A histogram of the bead score
#' should show a clear, small peak that represents the beads.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' sce <- initialBead(sce)
#' head(sce$scores)
#' head(sce$initial)
#'
#' @export
initialBead <- function(x) {
    
    unclassified.ind <- which(x$label == "cell")
    bead_channels <- grep("Bead", colnames(x$tech))
    
    if (min(as.matrix(x$tech[, bead_channels]), na.rm = TRUE) < 0) {
        stop("Bead data should all be non-negative")
    }
    
    if (length(bead_channels) > 1) {
        x$scores[, "beadScore"] <- rowSums(as.matrix(x$tech[, bead_channels]), 
                                           na.rm = FALSE)
    } else if (length(bead_channels) == 1){
        x$scores[, "beadScore"] <- x$tech[, bead_channels]
    } else {
        stop("No bead channels in data")
    }
    
    g <- initialGuess(x$scores$beadScore[unclassified.ind], middleGroup = 0)
    x$initial[unclassified.ind, "beadInitial"] <- g$label
    
    x
}
