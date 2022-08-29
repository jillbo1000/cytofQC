#' Preliminary doublet classification
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}}. 
#' @param score A value of 1, 2, or 3 that specifies the doublet score that 
#' should be calculated. See details for information on the doublet score. 
#' @param standardize A value of TRUE will use compute the doublet score 
#' using standardized data. The raw data will be used to compute the score 
#' if FALSE. It is highly recommended that the data are standardized prior
#' to computing the score because the variables are on different scales.
#'
#' @return A \code{SingleCellExperiment} that contains the doublet score and 
#' the doublet designation for each event. This information is stored in the 
#' \code{score} and \code{initial} objects in the colData for the 
#' \code{SingleCellExperiment}. 
#'
#' @details
#' The beads are typically the first cell classification that is done
#' because their identification is straightforward. Debris is typically
#' classified after the beads. This is because classifying debris
#' is more straightforward than doublets and labeling them before the
#' doublets aids in doublet classification. 
#'
#' Different event types are labeled iteratively so the \code{labels}
#' vector in the colData will contain all of the labels and 
#' probabilities computed up to this point. Only events that 
#' have a "cell" label can be assigned an initial event classification of
#' "doublet". This function computes a score that assesses how much an event
#' looks like a doublet and then fits a mixture model to assign each event 
#' a class of 1 for doublet, -1 for an event that is not a doublet, or 0 
#' for undetermined or previously assigned to a different event type. 
#' The score is recorded in the \code{score} object in the colData and 
#' the initial classification is recorded in the \code{initial} part of 
#' the colData. 
#' 
#' Several options are available for computing the doublet score. The
#' following list shows the doublet score calculations. Each one can
#' be selected by its number on the following list:
#' \enumerate{
#'   \item{DNA1 + DNA2 + Residual + Event_length - Offset - 0.5(Width)}
#'   \item{DNA1 + DNA2 + Residual + Event_length - Offset - 0.5(Width) + abs(Center)}
#'   \item{0.3 * (DNA1 + DNA2 + Event_length) + Residual + Center + (max(Offset))}
#' }
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' sce <- initialBead(sce)
#' sce <- initialDebris(sce)
#' sce <- initialDoublet(sce)
#' head(sce$scores)
#' head(sce$initial)
#'
#' @export
initialDoublet <- function(x, score = 1, standardize = TRUE) {
    
    if (standardize) {
        xs <- scale(x$tech)
    } else {
        xs <- x$tech
    }
    
    unclassified.ind <- which(x$label == "cell")
    
    if (score == 1) {
        x$scores[, "doubletScore"] <- xs[, "DNA1"] + xs[, "DNA2"] + 
            xs[, "Residual"] + xs[, "Event_length"] - xs[, "Offset"] - 
            0.5 * xs[, "Width"]
    } else if (score == 2) {
        x$scores[, "doubletScore"] <- xs[, "DNA1"] + xs[, "DNA2"] + 
            xs[, "Residual"] + xs[, "Event_length"] - xs[, "Offset"] - 
            0.5 * xs[, "Width"] + abs(xs[, "Center"])
    } else if (score == 3) {
        x$scores[, "doubletScore"] <- 0.3 * (xs[, "DNA1"] + xs[, "DNA2"] +
                                                 xs[, "Event_length"]) + 
            xs[, "Residual"] + xs[, "Center"] + 
            (max(xs[, "Offset"]) - xs[, "Offset"])
    } else {
        stop("Invalid score selection")
    }
    
    g <- initialGuess(x$scores$doubletScore[unclassified.ind], middleGroup = 1)
    x$initial$doubletInitial[unclassified.ind] <- g$label
    
    x
}