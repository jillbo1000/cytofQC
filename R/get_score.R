#' @name get_score
#' @rdname get_score
#' 
#' @title Returns a specified object from the cytofQC SingleCellExperiment
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}} 
#' with the scores and initial columns filled out for the event type of 
#' interest. 
#' @param type Identifies the type of object to be returned. for 
#' \code{get_score} and \code{get_prob}, type can be 'bead', 'dead', 'debris', 
#' or 'doublet' and it will return the numeric vector for the score or 
#' probability for the specified event type. If it is not specified, 
#' the \code{DataFrame} containing all of the scores or all of the 
#' probability will be returned. This argument does nothing for 
#' \code{get_labels} and \code{get_tech}. 
#'
#' @return For \code{get_prob} and \code{get_score}, a numeric vector 
#' containing the score or probabilities for the event type is returned. 
#' If \code{type} is not specified, a \code{DataFrame} containing all 
#' of the probabilities or scores is returned. 
#' 
#' For \code{get_labels}, a character vector containing the label
#' for each event is returned.
#' 
#' For \code{get_tech}, a \code{DataFrame} containing the technical variables
#' used to determine the label of each event. The bead, DNA, and viability
#' variables have an arcsinh transform, Event_length is unchanged, and 
#' the Gaussian parameters have a log transform using \code{log1p}.
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- initialBead(sce)
#' sce <- svmLabel(sce, type = "bead", loss = "auc")
#' cytofHist(get_score(sce, type = 'bead'), get_labels(sce), title = "Bead score")
NULL

#' @rdname get_score
#' @export
get_score <- function(x, type = NULL) {
    
    if (!is.null(type)) type <- tolower(type)

    if (!is.null(type)) {
        if(!(type %in% c("bead", "dead", "debris", "doublet"))) {
            stop("Invalid type. Must be NULL, 'bead', 'dead', 'debris', or 'doublet'.")
        }
    }
    
    if (is.null(type)) {
        x$scores
    } else {
        name <- paste0(type, "Score")
        x$scores[, name]
    }
}    

#' @rdname get_score
#' @export
get_prob <- function(x, type = NULL) {

    if (!is.null(type)) type <- tolower(type)
    
    if (!is.null(type)) {
        if(!(type %in% c("bead", "dead", "debris", "doublet"))) {
            stop("Invalid type. Must be NULL, 'bead', 'dead', 'debris', or 'doublet'.")
        }
    }
    
    if (is.null(type)) {
        x$probs
    } else {
        x$probs[, type]
    }
}

#' @rdname get_score
#' @export
get_labels <- function(x) {
    x$label
}

#' @rdname get_score
#' @export
get_tech <- function(x) {
    x$tech
}
