#' Preliminary doublet classification
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}}. 
#' @param score A value of "simple", "complex", "1" that specifies the doublet 
#' score that should be calculated. See details for information on the doublet 
#' score. 
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
#' \itemize{
#'   \item{simple: DNA + Residual + 0.5 * Event_length}
#'   \item{complex: DNA + Residual + 0.5 * Event_length - 
#'   0.5 * scale(abs(center) + abs(Width) + abs(Offset))}
#'   \item{1: 2 * DNA + Residual + Event_length - Offset - 0.5 * Width}
#' }
#' Score "1" is from the original release of this package and it works well for
#' doublet determination. The original scores "2" and "3" did not work well on
#' further examination and are removed from this version of the package. 
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' sce <- initialBead(sce)
#' sce <- initialDebris(sce)
#' sce <- initialDoublet(sce)
#' head(scores(sce))
#' head(initial(sce))
#'
#' @export
initialDoublet <- function(x, score = c("simple", "complex", "1"), 
                           standardize = TRUE) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    score <- match.arg(as.character(score), c("simple", "complex", "1"))

    dna_channels <- grep("DNA", colnames(x$tech))
    
    if (length(dna_channels) > 1) {
        DNA <- rowSums(as.matrix(x$tech[, dna_channels]), 
                       na.rm = FALSE)
    } else if (length(dna_channels) == 1){
        DNA <- x$tech[, dna_channels]
    } else {
        stop("No DNA channels in data")
    }
    
    xs <- x$tech
    xs$DNA <- DNA
    
    if (standardize) {
        xs <- scale(xs)
    } 
    
    unclassified.ind <- which(x$label == "cell")
    
    if (score == "1") {
        x$scores[, "doubletScore"] <- 2 * xs[, "DNA"] + xs[, "Residual"] + 
            xs[, "Event_length"] - xs[, "Offset"] - 0.5 * xs[, "Width"]
    } else if (score == "simple") {
        x$scores[, "doubletScore"] <- xs[, "DNA"] + 2 * xs[, "Residual"] + 
            0.5 * xs[, "Event_length"]
    } else if (score == "complex") {
        x$scores[, "doubletScore"] <- xs[, "DNA"] + xs[, "Residual"] + 
            0.5 * xs[, "Event_length"] + 
            0.5 * scale(abs(xs[, "Center"]) + abs(xs[, "Width"]) + abs(xs[, "Offset"]))
    }
    
    g_doublet <- initialGuess(x$scores$doubletScore[unclassified.ind], 
                              middleGroup = 1)$label
    g_normal <- initialGuess(x$scores$doubletScore[unclassified.ind], 
                             middleGroup = 0)$label
    g <- ifelse(g_doublet == 0, g_normal, g_doublet)
    x$initial$doubletInitial[unclassified.ind] <- g
    
    x
}