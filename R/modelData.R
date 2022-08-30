#' Returns indices for data to be used to create the final classification model
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}} 
#' with the scores and initial columns filled out for the event type of 
#' interest. 
#' @param type Identifies the type of label that is being modeled. Must
#' be 'bead', 'doublet', 'debris', or 'dead'.
#' @param n number of indices to return.
#'
#' @return An integer vector that contains the indices of the events that
#' should be included in the creation of the final classification model for
#' the event type of interest (bead, debris, doublet, dead).
#'
#' @details
#' The indices that are returned by \code{modelData} are be used to
#' create a model that can be used to classify the observations with
#' regard to the parameter of interest (bead, doublet, debris, dead).
#' It is used as part of \code{gbmLabel}, \code{rfLabel}, 
#' \code{svmLable}, and \code{labelQC}. The function \code{mocelData} 
#' uses the score and the function \code{initialGuess} to randomly select
#' a set of data points that we are confident are of the event type and 
#' not of the selected event type that can be used to train the data. Only
#' points that are labeled as \code{-1} and \code{1} are considered for the 
#' training dataset. The selected dataset is balance with a fairly equal 
#' number of points from each group. 
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- initialBead(sce)
#' train <- modelData(sce, type = "bead", n = 4000)
#'
#' @export
modelData <- function(x, type = c("bead", "doublet", "debris", "dead"), 
                      n = 4000) {
    
    type <- tolower(type)
    if (!(type %in% c("bead", "doublet", "debris", "dead"))) {
        stop("type must be either 'bead', 'doublet', 'debris', or 'dead'.")
    }
    
    if (length(type) != 1) {
        stop("Only one type can be selected.")
    }
    
    poss.ind <- seq_along(x$label)
    poss.ind <- poss.ind[x$initial[, grep(type, colnames(x$initial))] != 0 & 
                             x$label == "cell"]
    if (length(poss.ind) < n * 2) {
        n <- 0.5 * length(poss.ind)
        warning("Fewer than n/2 points in dataset. ", n, 
                " points used in training set.")
    }
    
    poss.wt <- ifelse(x$initial[, grep(type, 
                                       colnames(x$initial))][poss.ind] == -1, 
                      (1000 / 
                           table(x$initial[, 
                                           grep(type, 
                                                colnames(x$initial))][poss.ind]))[1], 
                      (1000 / 
                           table(x$initial[, 
                                           grep(type, 
                                                colnames(x$initial))][poss.ind]))[2])
    
    sample(poss.ind, n, prob = poss.wt)
    
}
