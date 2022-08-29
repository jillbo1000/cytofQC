#' Returns the final label assignments the specified parameters
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}} 
#' with the scores and initial columns filled out for the event type of 
#' interest. 
#' @param model Type of model to use to do the labeling. Options are 
#' "svm" for a support vector machine, "gbm" for a gradient boosting
#' machine, or "rf" for a random forest.
#' @param types Types of to model. Options are "bead", "doublet", 
#' "debris", and "dead".
#' @param nTrain The (maximum) number of data points to use when training a
#'   model to predict event types.
#' @param loss Specifies the type of loss used to tune the GBM. Can be either
#' "auc" for the area under the curve or "class" for classification error. 
#' This argument is ignored if random forest is used as the model. 
#'
#' @return A \code{SingleCellExperiment} data.frame is returned with the 
#' labels for the parameters of listed in \code{types} (bead, doublet, debris, 
#' or dead) added to the \code{label} variable and the probabilities for 
#' each of the columns pertaining to the parameters listed in \code{probs}.
#'
#' @details
#' \code{labelQC} uses a support vector machine, gradient boosting machine, 
#' or a random forest to compute the final labels
#' for the specified parameter types (bead, doublet, debris, or dead). 
#' The predicted probabilities for all of the observations are stored in 
#' the variable associated with that type for further analysis. Thus, it 
#' is possible to have a probability greater than 0.5 for 'debris' but still 
#' have a label of 'bead' if an observation was classified as a bead prior to
#' classifying the debris.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- labelQC(sce)
#'
#' @export
labelQC <- function(x, model = "svm", 
                    types = c("bead", "doublet", "debris", "dead"), 
                    nTrain = 4000, loss = "auc") {
    
    types <- tolower(types)
    if (length(setdiff(types, c("bead", "doublet", "debris", "dead")))) {
        stop("types must be either 'bead', 'doublet', 'debris', or 'dead'.")
    }
    
    xs <- scale(x$tech)
    
    loss <- tolower(loss)
    if (loss != "auc" & loss != "class") {
        warning("Invalid loss specified. AUC used to tune model.")
        loss <- "auc"
    }
    
    model <- tolower(model)
    
    if (model == "svm") {
        if ("bead" %in% types) {
            try(x <- initialBead(x), silent = TRUE)
            try(x <- svmLabel(x, type = "bead", loss = loss), silent = TRUE)
            x$scores$beadScore[is.na(x$scores$beadScore)] <- 0
            x$probs$bead[is.na(x$probs$bead)] <- 0
        } 
        
        if ("debris" %in% types) {
            try(x <- initialDebris(x), silent = TRUE)
            try(x <- svmLabel(x, type = "debris", loss = loss), silent = TRUE)
            x$scores$debrisScore[is.na(x$scores$debrisScore)] <- 0
            x$probs$debris[is.na(x$probs$debris)] <- 0
        }
        
        if ("doublet" %in% types) {
            try(x <- initialDoublet(x), silent = TRUE)
            try(x <- svmLabel(x, type = "doublet", loss = loss), silent = TRUE)
            x$scores$doubletScore[is.na(x$scores$doubletScore)] <- 0
            x$probs$doublet[is.na(x$probs$doublet)] <- 0
        }
        
        if ("dead" %in% types) {
            try(x <- initialDead(x), silent = TRUE)
            try(x <- svmLabel(x, type = "dead", loss = loss), silent = TRUE)
            x$scores$deadScore[is.na(x$scores$deadScore)] <- 0
            x$probs$dead[is.na(x$probs$dead)] <- 0
        }
        
    } else if (model == "gbm") {
        if ("bead" %in% types) {
            try(x <- initialBead(x), silent = TRUE)
            try(x <- gbmLabel(x, type = "bead", loss = loss), silent = TRUE)
            x$scores$beadScore[is.na(x$scores$beadScore)] <- 0
            x$probs$bead[is.na(x$probs$bead)] <- 0
        } 
        
        if ("debris" %in% types) {
            try(x <- initialDebris(x), silent = TRUE)
            try(x <- gbmLabel(x, type = "debris", loss = loss), silent = TRUE)
            x$scores$debrisScore[is.na(x$scores$debrisScore)] <- 0
            x$probs$debris[is.na(x$probs$debris)] <- 0
        }
        
        if ("doublet" %in% types) {
            try(x <- initialDoublet(x), silent = TRUE)
            try(x <- gbmLabel(x, type = "doublet", loss = loss), silent = TRUE)
            x$scores$doubletScore[is.na(x$scores$doubletScore)] <- 0
            x$probs$doublet[is.na(x$probs$doublet)] <- 0
        }
        
        if ("dead" %in% types) {
            try(x <- initialDead(x), silent = TRUE)
            try(x <- gbmLabel(x, type = "dead", loss = loss), silent = TRUE)
            x$scores$deadScore[is.na(x$scores$deadScore)] <- 0
            x$probs$dead[is.na(x$probs$dead)] <- 0
        }
        
    } else if (model == "rf") {
        if ("bead" %in% types) {
            try(x <- initialBead(x), silent = TRUE)
            try(x <- rfLabel(x, type = "bead"), silent = TRUE)
            x$scores$beadScore[is.na(x$scores$beadScore)] <- 0
            x$probs$bead[is.na(x$probs$bead)] <- 0
        } 
        
        if ("debris" %in% types) {
            try(x <- initialDebris(x), silent = TRUE)
            try(x <- rfLabel(x, type = "debris"), silent = TRUE)
            x$scores$debrisScore[is.na(x$scores$debrisScore)] <- 0
            x$probs$debris[is.na(x$probs$debris)] <- 0
        }
        
        if ("doublet" %in% types) {
            try(x <- initialDoublet(x), silent = TRUE)
            try(x <- rfLabel(x, type = "doublet"), silent = TRUE)
            x$scores$doubletScore[is.na(x$scores$doubletScore)] <- 0
            x$probs$doublet[is.na(x$probs$doublet)] <- 0
        }
        
        if ("dead" %in% types) {
            try(x <- initialDead(x), silent = TRUE)
            try(x <- rfLabel(x, type = "dead"), silent = TRUE)
            x$scores$deadScore[is.na(x$scores$deadScore)] <- 0
            x$probs$dead[is.na(x$probs$dead)] <- 0
        }
    } else {
        stop("Invalid model type")
    }
    
    x  
}
