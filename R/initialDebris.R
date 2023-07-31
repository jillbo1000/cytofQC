#' Preliminary debris classification
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}}. 
#' @param score A value of "simple", "complex", or "1" that specifies the debris 
#' score that should be calculated. Note that the score "1" is from the first
#' release of this package. See details for information on the debris score. 
#' @param standardize A value of TRUE will compute the debris score 
#' using standardized data. The raw data will be used to compute the score 
#' if FALSE. It is highly recommended that the data are standardized prior
#' to computing the score because the variables are on different scales.
#'
#' @return A \code{SingleCellExperiment} that contains the debris score and the 
#' debris designation for each event. This information is stored in the 
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
#' "debris". This function computes a score that assesses how much an event
#' looks like debris and then fits a mixture model to assign each event 
#' a class of 1 for debris, -1 for an event that is not debris, or 0 
#' for undetermined or previously assigned to a different event type. 
#' The score is recorded in the \code{score} object in the colData and 
#' the initial classification is recorded in the \code{initial} part of 
#' the colData. 
#' 
#' Several options are available for computing the debris score. The
#' following list shows the debris score calculations. Each one can
#' be selected by its number on the following list:
#' \itemize{
#'   \item{simple: -1 * (2 * DNA + Event_length)}
#'   \item{complex: -1 * (2 * DNA + Event_length - Offset - Width)}
#'   \item{1: 1 - (2 * DNA + Event_length - Center - Width + Offset)}
#' }
#' Note that there were three scores released with the original version of this
#' package that were labeled "1", "2", and "3". Score "1" is the same score and 
#' works fairly well. Scores "2" and "3" were found to be poor scores based
#' on further analysis so they are not included in this version of the package.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' sce <- initialBead(sce)
#' sce <- initialDebris(sce)
#' head(scores(sce))
#' head(initial(sce))
#'
#' @export
initialDebris <- function(x, score = c("simple", "complex", 1), standardize = TRUE) {
    
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
        x$scores[, "debrisScore"] <- 1 - (2 * xs[, "DNA"] +  
                                              xs[, "Event_length"] - 
                                              xs[, "Center"] - 
                                              xs[, "Width"] + 
                                              xs[, "Offset"])
    } else if (score == "simple") {
        x$scores[, "debrisScore"] <- -1 * (2 * xs[, "DNA"] + 
                                               xs[, "Event_length"])
        
    } else if (score == "complex") {
        x$scores[, "debrisScore"] <- -1 * (2 * xs[, "DNA"] + 
                                               xs[, "Event_length"]) -
            scale(xs[, "Offset"] - xs[, "Width"])[, 1]
    } 
    
    g_debris <- initialGuess(x$scores$debrisScore[unclassified.ind], middleGroup = 1)$label
    g_normal <- initialGuess(x$scores$debrisScore[unclassified.ind], middleGroup = 0)$label
    g <- ifelse(g_debris == 0, g_normal, g_debris)
    x$initial$debrisInitial[unclassified.ind] <- g
    
    x
}
