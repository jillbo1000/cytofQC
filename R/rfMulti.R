#' Returns the final label assignments for a parameter using a multiclass
#' random forest
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
#' \code{rfLabel} uses a random forest to compute the final labels
#' for the specified parameter type (bead, doublet, debris, or dead). The
#' model is computed using only the data specified in the index argument. The
#' random forest model is computed using the default settings in
#' \code{randomForest} using the observations specified in \code{index}.
#' The predicted values are computed
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
#' rfBeads <- rfLabel(x, labels, type = 'bead', init = beads$init, index = index)
#'
#' @export
rfLabel <- function(x, labels, type, init, index, standardize = TRUE) {
  
  type <- tolower(type)
  if (!(type %in% c("bead", "doublet", "debris", "dead"))) {
    stop("type must be either 'bead', 'doublet', 'debris', or 'dead'.")
  }
  
  Time <- subset(x, select = c(get("Time")))
  x <- subset(x, select = -c(get("Time")))
  
  if (standardize) {
    x <- as.data.frame(scale(x))
  }
  
  rffit <- randomForest::randomForest(x = x[index, ], y = factor(init[index]))
  pred <- stats::predict(rffit, x, type = 'prob')[, 2]
  
  labs <- labels
  labs[, type] <- pred
  labs$label[labs$label == "cell"] <- ifelse(round(pred[labs$label == "cell"]),
                                             type, 'cell')
  
  labs
  
}
