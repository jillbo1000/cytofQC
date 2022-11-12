#' @name scores
#' @rdname scores
#' 
#' @title Returns a specified object from the cytofQC SingleCellExperiment
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}} 
#' with the scores and initial columns filled out for the event type of 
#' interest. 
#' @param type Identifies the type of objects to be returned. For 
#' \code{scores} and \code{probs}, type can be one or more of 'bead', 
#' 'dead', 'debris', or 'doublet'. For \code{tech} it can be any of the QC 
#' variables. It will return the numeric vector or \code{DataFrame} 
#' for the score(s) or probability for the specified event type(s). If the 
#' event types are not specified, the \code{DataFrame} containing all of the 
#' scores or all of the probability will be returned. This argument does 
#' nothing for \code{label}. 
#'
#' @return For \code{probs}, \code{scores}, and \code{tech}, a numeric vector 
#' or \code{DataFrame} with the information for the event type(s) is returned. 
#' 
#' For \code{label}, a character vector containing the label for each event is 
#' returned.
#' 
#' For \code{tech}, a \code{DataFrame} containing the technical variables
#' used to determine the label of each event. The bead, DNA, and viability
#' variables have an arcsinh transform, Event_length is unchanged, and 
#' the Gaussian parameters have a log transform using \code{log1p}.
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- labelQC(sce)
#' table(label(sce))
#' cytofHist(scores(sce, type = 'bead'), label(sce), title = "Bead score")
NULL

#' @rdname scores
#' @export
scores <- function(x, type = c("all", "bead", "debris", "doublet", "dead")) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    type <- match.arg(tolower(type), choices = c("all", "bead", "debris", 
                                                 "doublet", "dead"), 
                      several.ok = TRUE)
    
    if ("all" %in% type) {
        x$scores
    } else {
        nm <- paste0(type, "Score")
        x$scores[, nm]
    }
}    

#' @rdname scores
#' @export
probs <- function(x, type = c("all", "bead", "debris", "doublet", "dead")) {

    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    type <- match.arg(tolower(type), choices = c("all", "bead", "debris", 
                                                 "doublet", "dead"), 
                      several.ok = TRUE)
    
    if ("all" %in% type) {
        x$probs
    } else {
        nm <- colnames(x$probs)[colnames(x$probs) %in% type]
        x$probs[, nm]
    }
}

#' @rdname scores
#' @export
label <- function(x) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    x$label
}

#' @rdname scores
#' @export
tech <- function(x, type = c("all", "Bead", "DNA", "Viability", 
                               "Event_length", "Center", "Offset", "Width", 
                               "Residual")) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    type <- match.arg(tolower(type), choices = c("all", "Bead", "DNA", 
                                                 "Viability", "Event_length", 
                                                 "Center", "Offset", "Width", 
                                                 "Residual"), 
                      several.ok = TRUE)
    
    if ("all" %in% type) {
        x$tech
    } else {
        nm1 <- tolower(type)
        nm2 <- tolower(colnames(x$tech))
        pattern <- paste(nm1, collapse="|")
        ind <- grep(pattern, nm2)
        x$tech[, ind]
    }
}

#' @rdname scores
#' @export
initial <- function(x, type = c("all", "bead", "debris", "doublet", "dead")) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    type <- match.arg(tolower(type), choices = c("all", "bead", "debris", 
                                                 "doublet", "dead"), 
                      several.ok = TRUE)
    
    if ("all" %in% type) {
        x$initial
    } else {
        nm <- paste0(type, "Initial")
        x$initial[, nm]
    }
    
}

