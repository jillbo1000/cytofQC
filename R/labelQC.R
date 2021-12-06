#' Returns the final label assignments the specified parameters
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}.
#' @param model Type of model to use to do the labeling. Options are 
#' "svm" for a support vector machine, "gbm" for a gradient boosting
#' machine, or "rf" for a random forest.
#' @param types Types of to model. Options are "bead", "doublet", 
#' "debris", and "dead". 
#' @param loss Specifies the type of loss used to tune the model. Can 
#' be either "auc" or "class". This argument is ignored if random
#' forest is used as the model. 
#'
#' @return A \code{label} data.frame is returned with the labels
#' for the parameters of listed in \code{types} (bead, doublet, debris, 
#' or dead) added to the \code{label} variable and the probabilities for 
#' each of the columns pertaining to the parameter listed in \code{types}.
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
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' svmLabels <- labelQC(x, model = "svm", types = c("bead", "doublet", "debris"))
#'
#' @export
labelQC <- function(x, model = "svm", types = c("bead", "doublet", "debris", "dead"), loss = "auc") {
  
  types <- tolower(types)
  if (length(setdiff(types, c("bead", "doublet", "debris", "dead")))) {
    stop("types must be either 'bead', 'doublet', 'debris', or 'dead'.")
  }
  
  labels <- qcDataFrame(x)
  
  Time <- x[, 1]
  x <- x[, -1]
  x <- as.data.frame(scale(x))
  
  loss <- tolower(loss)
  if (loss != "auc" & loss != "class") {
    warning("Invalid loss specified. AUC used to tune model.")
    loss <- "auc"
  }
  
  mnn <- NN(x, n = 20, type = "mnn")
  
  model <- tolower(model)
  
  if (model == "svm") {
    if ("bead" %in% types) {
      beads <- initialBead(x, labels)
      sure <- surematch(mnn, beads$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = beads$init, n = 4000)
      labels <- svmLabel(x, labels, type = "bead", init = beads$init, index = ind, loss = loss)
    } 
    
    if ("doublet" %in% types) {
      doublets <- initialDoublet(x, labels, score = 1)
      sure <- surematch(mnn, doublets$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = doublets$init, n = 4000)
      labels <- svmLabel(x, labels, type = "doublet", init = doublets$init, index = ind, loss = loss)
    }
    
    if ("debris" %in% types) {
      debris <- initialDebris(x, labels, score = 1)
      sure <- surematch(mnn, debris$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = debris$init, n = 4000)
      labels <- svmLabel(x, labels, type = "debris", init = debris$init, index = ind, loss = loss)
    }
    
    if ("dead" %in% types) {
      dead <- initialDead(x, labels, dna = TRUE)
      sure <- surematch(mnn, dead$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = dead$init, n = 4000)
      labels <- svmLabel(x, labels, type = "dead", init = dead$init, index = ind, loss = loss)
    }	  
  } else if (model == "gbm") {
    if ("bead" %in% types) {
      beads <- initialBead(x, labels)
      sure <- surematch(mnn, beads$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = beads$init, n = 4000)
      labels <- gbmLabel(x, labels, type = "bead", init = beads$init, index = ind, loss = loss)
    } 
    
    if ("doublet" %in% types) {
      doublets <- initialDoublet(x, labels, score = 1)
      sure <- surematch(mnn, doublets$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = doublets$init, n = 4000)
      labels <- gbmLabel(x, labels, type = "doublet", init = doublets$init, index = ind, loss = loss)
    }
    
    if ("debris" %in% types) {
      debris <- initialDebris(x, labels, score = 1)
      sure <- surematch(mnn, debris$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = debris$init, n = 4000)
      labels <- gbmLabel(x, labels, type = "debris", init = debris$init, index = ind, loss = loss)
    }
    
    if ("dead" %in% types) {
      dead <- initialDead(x, labels, dna = TRUE)
      sure <- surematch(mnn, dead$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = dead$init, n = 4000)
      labels <- gbmLabel(x, labels, type = "dead", init = dead$init, index = ind, loss = loss)
    }	  
  } else if (model == "rf") {
    if ("bead" %in% types) {
      beads <- initialBead(x, labels)
      sure <- surematch(mnn, beads$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = beads$init, n = 4000)
      labels <- rfLabel(x, labels, type = "bead", init = beads$init, index = ind)
    } 
    
    if ("doublet" %in% types) {
      doublets <- initialDoublet(x, labels, score = 1)
      sure <- surematch(mnn, doublets$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = doublets$init, n = 4000)
      labels <- rfLabel(x, labels, type = "doublet", init = doublets$init, index = ind)
    }
    
    if ("debris" %in% types) {
      debris <- initialDebris(x, labels, score = 1)
      sure <- surematch(mnn, debris$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = debris$init, n = 4000)
      labels <- rfLabel(x, labels, type = "debris", init = debris$init, index = ind)
    }
    
    if ("dead" %in% types) {
      dead <- initialDead(x, labels, dna = TRUE)
      sure <- surematch(mnn, dead$init, threshold = 0)
      ind <- modelData(labels, subset = sure$match, init = dead$init, n = 4000)
      labels <- rfLabel(x, labels, type = "dead", init = dead$init, index = ind)
    }	    
  } else {
    stop("Invalid model type")
  }
  
  labels  
}
