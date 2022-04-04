#' Returns the final label assignments the specified parameters
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}.
#' @param model Type of model to use to do the labeling. Options are 
#' "svm" for a support vector machine, "gbm" for a gradient boosting
#' machine, or "rf" for a random forest.
#' @param types Types of to model. Options are "bead", "doublet", 
#' "debris", and "dead".
#' @param nn Specifies if k-nearest neighbors ("knn") or mutual nearest
#' neighbors ("mnn") for initial sample selection. 
#' @param nnNum The number of nearest neighbors to compute in the 
#' nearest neighbors matrix. 
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
labelQC <- function(x, model = "svm", types = c("bead", "doublet", "debris", "dead"), 
                    nTrain = 4000, loss = "auc") {
    
    types <- tolower(types)
    if (length(setdiff(types, c("bead", "doublet", "debris", "dead")))) {
        stop("types must be either 'bead', 'doublet', 'debris', or 'dead'.")
    }
    
    labels <- qcDataFrame(x)
    17
    Time <- subset(x, select = c(get("Time")))
    x <- subset(x, select = -c(get("Time")))
    x <- as.data.frame(scale(x))
    x <- cbind(Time, x)
    
    tTest <- nTrain / 2
    
    loss <- tolower(loss)
    if (loss != "auc" & loss != "class") {
        warning("Invalid loss specified. AUC used to tune model.")
        loss <- "auc"
    }
    
    model <- tolower(model)
    
    if (model == "svm") {
        if ("bead" %in% types) {
            beads <- initialBead(x, labels)
            if (min(sum(beads$init == -1), sum(beads$init == 1)) < tTest) {
                types <- types[types != "bead"]
                warning("Not enough beads or non-beads in dataset to train a model with nTrain value specified. Bead data not fitted.")
            } else {
                ind <- modelData(labels, init = beads$init , n = nTrain)
                labels <- svmLabel(x, labels, type = "bead", init = beads$init, index = ind, loss = loss)
            }
        } 
        
        if ("debris" %in% types) {
            debris <- initialDebris(x, labels, score = 1)
            if (min(sum(debris$init==-1), sum(debris$init==1)) < tTest) {
                types <- types[types != "debris"]
                warning("Not enough debris or non-debris observations in dataset to train a model with nTrain value specified. Debris not fitted.")
            } else {
                ind <- modelData(labels, init = debris$init, n = nTrain)
                labels <- svmLabel(x, labels, type = "debris", init = debris$init, index = ind, loss = loss)
            }
        }
        
        if ("doublet" %in% types) {
            doublets <- initialDoublet(x, labels, score = 1)
            if (min(sum(doublets$init == -1), sum(doublets$init == 1)) < tTest) {
                types <- types[types != "doublet"]
                warning("Not enough doublets or non-doublet observations in dataset to train a model with nTrain value specified. Doublets not fitted.")
            } else {
                ind <- modelData(labels, init = doublets$init, n = nTrain)
                labels <- svmLabel(x, labels, type = "doublet", init = doublets$init, index = ind, loss = loss)
            }
        }
        
        if ("dead" %in% types) {
            dead <- initialDead(x, labels, dna = TRUE)
            if (min(sum(dead$init == -1), sum(dead$init == 1)) < tTest) {
                types <- types[types != "dead"]
                warning("Not enough dead or non-dead observations in dataset to train a model with nTrain value specified. Live/dead not fitted.")
            } else {
                ind <- modelData(labels, init = dead$init, n = nTrain)
                labels <- svmLabel(x, labels, type = "dead", init = dead$init, index = ind, loss = loss)
            }
        }	  
    } else if (model == "gbm") {
        if ("bead" %in% types) {
            beads <- initialBead(x, labels)
            if (min(sum(beads$init == -1), sum(beads$init == 1)) < tTest) {
                types <- types[types != "bead"]
                warning("Not enough beads or non-beads in dataset to train a model with nTrain value specified. Bead data not fitted.")
            } else {
                ind <- modelData(labels, init = beads$init, n = nTrain)
                labels <- gbmLabel(x, labels, type = "bead", init = beads$init, index = ind, loss = loss)
            }
        } 
        
        if ("debris" %in% types) {
            debris <- initialDebris(x, labels, score = 1)
            if (min(sum(debris$init == -1), sum(debris$init == 1)) < tTest) {
                types <- types[types != "debris"]
                warning("Not enough debris or non-debris observations in dataset to train a model with nTrain value specified. Debris not fitted.")
            } else {
                ind <- modelData(labels, init = debris$init, n = nTrain)
                labels <- gbmLabel(x, labels, type = "debris", init = debris$init, index = ind, loss = loss)
            }
        }
        
        if ("doublet" %in% types) {
            doublets <- initialDoublet(x, labels, score = 1)
            if (min(sum(doublets$init == -1), sum(doublets$init == 1)) < tTest) {
                types <- types[types != "doublet"]
                warning("Not enough doublets or non-doublet observations in dataset to train a model with nTrain value specified. Doublets not fitted.")
            } else {
                ind <- modelData(labels, init = doublets$init, n = nTrain)
                labels <- gbmLabel(x, labels, type = "doublet", init = doublets$init, index = ind, loss = loss)
            }
        }
        
        if ("dead" %in% types) {
            dead <- initialDead(x, labels, dna = TRUE)
            if (min(sum(dead$init == -1), sum(dead$init == 1)) < tTest) {
                types <- types[types != "dead"]
                warning("Not enough dead or non-dead observations in dataset to train a model with nTrain value specified. Live/dead not fitted.")
            } else {
                ind <- modelData(labels, init = dead$init, n = nTrain)
                labels <- gbmLabel(x, labels, type = "dead", init = dead$init, index = ind, loss = loss)
            }
        }	  
    } else if (model == "rf") {
        if ("bead" %in% types) {
            beads <- initialBead(x, labels)
            if (min(sum(beads$init == -1), sum(beads$init == 1)) < tTest) {
                types <- types[types != "bead"]
                warning("Not enough beads or non-beads in dataset to train a model with nTrain value specified. Bead data not fitted.")
            } else {
                ind <- modelData(labels, init = beads$init, n = nTrain)
                labels <- rfLabel(x, labels, type = "bead", init = beads$init, index = ind)
            }
        } 
        
        if ("debris" %in% types) {
            debris <- initialDebris(x, labels, score = 1)
            if (min(sum(debris$init == -1), sum(debris$init == 1)) < tTest) {
                types <- types[types != "debris"]
                warning("Not enough debris or non-debris observations in dataset to train a model with nTrain value specified. Debris not fitted.")
            } else {
                ind <- modelData(labels, init = debris$init, n = nTrain)
                labels <- rfLabel(x, labels, type = "debris", init = debris$init, index = ind)
            }
        }
        
        if ("doublet" %in% types) {
            doublets <- initialDoublet(x, labels, score = 1)
            if (min(sum(doublets$init == -1), sum(doublets$init == 1)) < tTest) {
                types <- types[types != "doublet"]
                warning("Not enough doublets or non-doublet observations in dataset to train a model with nTrain value specified. Doublets not fitted.")
            } else {
                ind <- modelData(labels, init = doublets$init, n = nTrain)
                labels <- rfLabel(x, labels, type = "doublet", init = doublets$init, index = ind)
            }
        }
        
        if ("dead" %in% types) {
            dead <- initialDead(x, labels, dna = TRUE)
            if (min(sum(dead$init == -1), sum(dead$init == 1)) < tTest) {
                types <- types[types != "dead"]
                warning("Not enough dead or non-dead observations in dataset to train a model with nTrain value specified. Live/dead not fitted.")
            } else {
                ind <- modelData(labels, init = dead$init, n = nTrain)
                labels <- rfLabel(x, labels, type = "dead", init = dead$init, index = ind)
            }
        }	  
    } else {
        stop("Invalid model type")
    }
    
    labels  
}
