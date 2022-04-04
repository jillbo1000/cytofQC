#' Preliminary viability classification.
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}. Note that the
#' data should not be standardized. All of the bead channels should be
#' non-negative.
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param dna If TRUE, the DNA will be used to determine viability.
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' this analysis.
#'
#' @return A \code{data.frame} that contains the viablity designation
#' for the observation. The first column is the \code{Time} variable
#' that helps ensure that all of the observations can be correctly matched
#' during all stages labeling. The last column is \code{init} which is
#' TRUE if the observation is determined to be dead and FALSE if the
#' observation is not dead during this initial classification.
#'
#' @details
#' The beads are typically the first cell classification that is done. The
#' different cell types are labeled iteratively so the \code{labels}
#' data.frame should contain all of the labels and probabilities computed
#' up to this point. Thus, if the dead cells were not the first cell type that
#' was classified, the \code{labels} data.frame must have the labels and
#' probabilities computed prior to classifying the dead cells. Also note that
#' this matrix will not be changed during this state. This function
#' creates the most preliminary bead labels and that information is needed
#' to obtain the final bead labels. The \code{labels} data.frame is needed
#' to determine which cells have and have not been labeled prior to this
#' step.
#'
#' The viability measures should classify into two fairly clear groups
#' where one is permeable cells and the other is non-permeability. DNA is
#' also often higher for cells that are permeable. The primary measure
#' for determining permeability is sum of the viability and the DNA measures.
#' The function \code{find_groups} is used to determine the groups. The members of the group with the
#' largest mean are classified as beads and the rest are classified as
#' not beads. The observations that are labeled 'GDPzero' are included in
#' the output data, but the are removed from the data for bead designation.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' tech <- dataPrep(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' labels <- qcDataFrame(tech)
#' initialDead(tech, labels)
#'
#' @export
initialDead <- function(x, labels, dna = TRUE, standardize = TRUE) {
    
    if (standardize) {
        xs <- scale(x[, -1])
    } else {
        xs <- x
    }
    
    unclassified.ind <- which(labels$label == "cell")

    viability <- xs[, grep("Viability", colnames(xs))]

    if (ncol(data.frame(viability)) > 1) {
        viability <- rowSums(viability, na.rm = TRUE)
    } 
    
    if (dna) {
        dnaData <- x[, grep("DNA", colnames(x))]
        if (ncol(dnaData) > 1) {
            n <- ncol(dnaData)
            dnaData <- rowSums(dnaData, na.rm = TRUE)
            dnaData <- dnaData / n
        }
        deadScore <- viability + 0.5 * dnaData
    } else {
        deadScore <- viability
    }
    
    g <- initialGuess(deadScore[unclassified.ind], middleGroup = 1)
    init <- rep(0, nrow(x))
    init[unclassified.ind] <- g$label
    
    data.frame(Time = x[, "Time"], deadScore = deadScore, init = init)
    
}
