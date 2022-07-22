#' Returns the final label assignments for a parameter using a semi-supervised
#' support vector machine
#' 
#' @param x A \code{SingleCellExperiment} created with \code{\link{readCytof}} 
#' with the scores and initial columns filled out for the event type of 
#' interest. 
#' @param type Identifies the type of label that is being modeled. Must
#' be 'bead', 'doublet', 'debris', or 'dead'.
#' @param loss Specifies the type of loss used to tune the SVM. Can be either
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
#' @details \code{s3vmLabel} uses a semi-supervised support vector machine to
#'   compute the final labels for the specified parameter type (bead, doublet,
#'   debris, or dead). The model is initially computed using only the data
#'   specified in the index argument. Events are iteratively added to this set
#'   when the updated SVM predicts a label with high confidence. Then predicted
#'   values are computed for all of the observations in \code{x}. If the
#'   predicted probability for the label type is greater than 0.5, the label is
#'   changed to the specified type. However, if an observation already has a
#'   label other than 'cell' in the \code{labels$label} variable, it will not be
#'   changed. The predicted probabilities for all of the observations is stored
#'   in the variable associated with that type for further analysis. Thus, it is
#'   possible to have a probability greater than 0.5 for 'debris' but still have
#'   a label of 'bead' if an observation was classified as a bead prior to
#'   classifying the debris.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- initialBead(sce)
#' sce <- svmLabel(sce, type = "bead", loss = "auc")
#'
#' @export
s3vmLabel <- function(x, type = c("bead", "doublet", "debris", "dead"), 
                      loss = "auc", n = 4000, standardize = TRUE) {
    
       
        type <- tolower(type)
        if (!(type %in% c("bead", "doublet", "debris", "dead"))) {
            stop("type must be either 'bead', 'doublet', 'debris', or 'dead'.")
        }
        
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
            warning(paste("Not enough ", type, " or non-", type, "to build model."))
            pred.pr <- rep(0, nrow(xs))
            
        } else {
            
            y.s3vm <- rep(NA, nrow(xs))
            y.s3vm[index] <- x$initial[index, grep(type, colnames(x$initial))]
            s3vmfit <- ssc::selfTraining(x = xs, y = factor(y.s3vm), 
                                         x.inst = TRUE, learner = e1071::svm,
                                         learner.pars = list(kernel ="radial", 
                                                             scale = FALSE, 
                                                             probability = TRUE),
                                         pred = function(m, x){
                                             attr(predict(m, x, probability = TRUE), 
                                                  "probabilities")
                                         })
            pred <- stats::predict(s3vmfit$model, xs, probability = TRUE)
            pred.pr <- attr(pred, "probabilities")[, colnames(attr(pred, "probabilities")) == "1"]
            
        }
        
        x$probs[, grep(type, colnames(x$initial))] <- pred.pr
        x$label[x$label == "cell"] <- ifelse(round(pred.pr[x$label == "cell"]), type, "cell")
        
        x
        
    }
    