#' Returns the final label assignments for a parameter using a support vector
#' machine
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}.
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param type Identifies the type of label that is being modeled. Must
#' be 'bead', 'doublet', 'debris', or 'dead'.
#' @param init A logical vector that contains the initial labeling for the
#' cells for the cell type of interest.
#' @param index A vector containing the indices of the data that should be
#' used to compute the model. These should be obtained from
#' \code{\link{modelData}}.
#' @param loss Specifies the type of loss used to tune the GBM. Can be either
#' "auc" or "class". 
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' this analysis.
#'
#' @return An updated \code{label} data.frame is returned with the labels
#' for the parameter of interest (bead, doublet, debris, or dead) added to
#' the \code{label} variable and the probabilities for the column
#' pertaining to the parameter filled in.
#'
#' @details
#' \code{svmLabel} uses a support vector machine to compute the final labels
#' for the specified parameter type (bead, doublet, debris, or dead). The
#' model is computed using only the data specified in the index argument. The svm
#' is tuned using \code{\link{EZtune}} and then predicted values are computed
#' for all of the observations in \code{x}. If the predicted probability for
#' the label type is greater than 0.5, the label is changed to the specified
#' type. However, if an observation already has a label other than 'cell'
#' in the \code{labels$label} variable, it will not be changed. The predicted
#' probabilities for all of the observations is stored in the variable
#' associated with that type for further analysis. Thus, it is possible
#' to have a probability greater than 0.5 for 'debris' but still have a
#' label of 'bead' if an observation was classified as a bead prior to
#' classifying the debris.
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' nn <- NN(x)
#' labels <- qcDataFrame(x)
#' beads <- initialBead(x, labels = labels)
#' mm <- mismatch(nn, init = beads$init, threshold = 1)
#' index <- modelData(labels, subset = mm$Mismatch, init = beads$init, n = 4000)
#' svmBeads <- svmLabel(x, labels, type = 'bead', init = beads$init, index = index)
#'
#' @export
svmLabel <- function(x, labels, type, init, index, loss = "auc", standardize = TRUE) {

  type <- tolower(type)
  if (!(type %in% c("bead", "doublet", "debris", "dead"))) {
    stop("type must be either 'bead', 'doublet', 'debris', or 'dead'.")
  }

  Time <- x[, 1]
  x <- x[, -1]

  if (standardize) {
    x <- as.data.frame(scale(x))
  }
  
  loss <- tolower(loss)
  if (loss != "auc" & loss != "class") {
    warning("Invalid loss specified. AUC used to tune model.")
    loss <- "auc"
  }
  
  svmTune <- EZtune::eztune(x = x[index, ], y = factor(init[index]),
                            method = "svm", fast = 0.5, loss = loss)
  svmfit <- e1071::svm(x = x[index, ], y = factor(init[index]),
                       cost = svmTune$cost, gamma = svmTune$gamma,
                       kernel = "radial", probability = TRUE)
  pred <- stats::predict(svmfit, x, probability = TRUE)
  pred.pr <- attr(pred, "probabilities")[, colnames(attr(pred, "probabilities")) == "TRUE"]

  labs <- labels
  labs[, type] <- pred.pr
  labs$label[labs$label == "cell"] <- ifelse(round(pred.pr[labs$label == "cell"]),
                                         type, 'cell')

  labs

}
