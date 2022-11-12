#' Preliminary viability classification
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}}. 
#' @param dna If TRUE, the DNA will be used to determine viability.
#' @param standardize A value of TRUE will use compute the doublet score 
#' using standardized data. The raw data will be used to compute the score 
#' if FALSE. It is highly recommended that the data are standardized prior
#' to computing the score because the variables are on different scales.
#'
#' @return A \code{SingleCellExperiment} that contains the permeability score 
#' and the permeability designation for each event. This information is stored 
#' in the \code{score} and \code{initial} objects in the colData for the 
#' \code{SingleCellExperiment}. 
#' 
#' @details
#' The beads are typically the first cell classification that is done
#' because their identification is straightforward. Debris is typically
#' classified after the beads and doublets classified after the debris. 
#' After the doublets are classified, the permeability is assessed. 
#' The permeability can be used to determine which cells are alive 
#' and which are dead. Dead cells are often gated out prior to gating
#' debris and doublets. However, cleaning with respect to permeability
#' is complicated in that it does not always make sense to clean "dead", 
#' or permeable, cells out of the data. Our default method chooses 
#' to classify them last because once an event is labeled "dead", it 
#' will not be assigned a different label. However, a user may wish to 
#' label the permeable, or dead cells right after cleaning the 
#' beads so that the "dead" label takes precedence over debris and
#' doublets. Note that "dead" is used in place of permeable for 
#' labeling in this package. 
#'
#' Different event types are labeled iteratively so the \code{labels}
#' vector in the colData will contain all of the labels and 
#' probabilities computed up to this point. Only events that  
#' have a "cell" label can be assigned an initial event classification of
#' "dead". This function computes a score that assesses how much an event
#' looks like a permeable cell and then fits a mixture model to assign 
#' each event a class of 1 for permeable, -1 for an event that is not 
#' permeable, or 0 for undetermined or previously assigned to a different 
#' event type. The score is recorded in the \code{score} object in the 
#' colData and the initial classification is recorded in the \code{initial} 
#' part of the colData. 
#' 
#' The viability measures should classify into two fairly clear groups
#' where one is permeable cells and the other is non-permeability. DNA is
#' also often higher for cells that are permeable. The primary measure
#' for determining permeability is sum of the viability measures, but 
#' the method allows for DNA content to be used as well.
#' The function \code{initialGuess} is used to determine the groups. The 
#' members of the group with the largest mean are classified as 'dead' 
#' and the rest are classified as not dead. 
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' sce <- initialBead(sce)
#' sce <- initialDebris(sce)
#' sce <- initialDoublet(sce)
#' sce <- initialDead(sce)
#' head(scores(sce))
#' head(initial(sce))
#'
#' @export
initialDead <- function(x, dna = FALSE, standardize = TRUE) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    if (standardize) {
        xs <- scale(x$tech)
    } else {
        xs <- x$tech
    }
    
    unclassified.ind <- which(x$label == "cell")
    
    viability <- xs[, grep("Viability", colnames(x$tech))]
    
    if (ncol(data.frame(viability)) > 1) {
        viability <- rowSums(viability, na.rm = TRUE)
    } 
    
    if (dna) {
        dnaData <- xs[, grep("DNA", colnames(x$tech))]
        if (ncol(dnaData) > 1) {
            n <- ncol(dnaData)
            dnaData <- rowSums(dnaData, na.rm = TRUE)
            dnaData <- dnaData / n
        }
        x$scores$deadScore <- viability + dnaData
    } else {
        x$scores$deadScore <- viability
    }
    
    g <- initialGuess(x$scores$deadScore[unclassified.ind], middleGroup = 1)
    x$initial[unclassified.ind, "deadInitial"] <- g$label
    
    x
}
