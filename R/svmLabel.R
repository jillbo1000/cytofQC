#' Returns the final label assignments for a parameter using a support vector
#' machine
#'
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}} 
#' with the scores and initial columns filled out for the event type of 
#' interest. 
#' @param type Identifies the type of label that is being modeled. Must
#' be 'bead', 'doublet', 'debris', or 'dead'.
#' @param loss Specifies the type of loss used to tune the GBM. Can be either
#' "auc" for the area under the curve or "class" for classification error. 
#' @param n number of observations in training dataset.
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' this analysis because the variables are on different scales.
#'
#' @return An updated \code{SingleCellExperiment} is returned with the labels
#' for the parameter of interest (bead, doublet, debris, or dead) added to
#' the \code{label} object of the \code{SingleCellExperiment} and the 
#' probabilities for the event type added to the \code{probs} object of the 
#' \code{SingleCellExperiment}.
#'
#' @details
#' \code{svmLabel} uses a support vector machine to compute the final labels
#' for the specified parameter type (bead, doublet, debris, or dead). This step
#' cannot be completed until the corresponding initialization function 
#' (\code{initialBead}, \code{initialDebris}, \code{initialDoublet}, or 
#' \code{initialDead}) is done on the \code{SingleCellExperiment} created by 
#' \code{readCytof}. 
#' The support vector machine is tuned using \code{\link[EZtune]{eztune}} and 
#' then predicted values are computed for all of the events in \code{x}. If 
#' the predicted probability for the label type is greater than 0.5, the label 
#' is changed to the specified type. However, if an observation already has a 
#' label other than 'cell' in the \code{label} variable, it will not be 
#' changed. The predicted probabilities for all of the observations are 
#' stored in the variable associated with that type in the \code{probs} 
#' object of \code{x} for further analysis. Thus, it is possible to have 
#' a probability greater than 0.5 for 'debris' but still have a label of 
#' 'bead' if an observation was classified as a bead prior to classifying 
#' the debris.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- initialBead(sce)
#' sce <- svmLabel(sce, type = "bead", loss = "auc")
#' table(label(sce))
#'
#' @export
svmLabel <- function(x, type = c("bead", "doublet", "debris", "dead"), 
                     loss = c("auc", "class"), n = 4000, standardize = TRUE) {
    
    if (!methods::is(x, "SingleCellExperiment")) {
        stop("x must be an object created with readCytof")
    }
    
    type <- match.arg(tolower(type), choices = c("bead", "debris", 
                                                 "doublet", "dead"))
    loss <- match.arg(tolower(loss), c("auc", "class"))
    
    if (standardize) {
        xs <- scale(x$tech)
    } else {
        xs <- x$tech
    }
    
    loss <- tolower(loss)
    if (loss != "auc" & loss != "class") {
        warning("Invalid loss specified. AUC used to tune model.")
        loss <- "auc"
    }
    
    index <- modelData(x, type = type)
    
    if (sum(x$initial[index, grep(type, colnames(x$initial))] == -1) < 100 | 
        sum(x$initial[index, grep(type, colnames(x$initial))] == 1) < 100) {
        warning("Not enough ", type, " or non-", type, "to build model.")
        pred <- rep(0, nrow(xs))
    } else {
        svmTune <- EZtune::eztune(x = xs[index, ], 
                                  y = factor(x$initial[index, 
                                                       grep(type, 
                                                            colnames(x$initial))]), 
                                  method = "svm",
                                  fast = 0.5, 
                                  loss = loss)
        
        svmfit <- e1071::svm(x = xs[index, ], 
                             y = factor(x$initial[index, 
                                                  grep(type, 
                                                       colnames(x$initial))]),
                             cost = svmTune$cost, 
                             gamma = svmTune$gamma,
                             kernel = "radial", 
                             probability = TRUE)
        
        pred.pr <- stats::predict(svmfit, xs, probability = TRUE)
        pred <- attr(pred.pr, 
                     "probabilities")[, 
                                      colnames(attr(pred.pr, 
                                                    "probabilities")) == "1"]
    }
    
    
    x$probs[, grep(type, colnames(x$initial))] <- pred
    x$label[x$label == "cell"] <- ifelse(round(pred[x$label == "cell"]), 
                                         type, "cell")
    
    x
}
