#' Returns indices for data to be used to create the final classification model.
#'
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param init A numeric vector that contains the initial labeling for the
#' observations for the parameter of interest.
#' @param subset A logical vector that indicates if an observation should
#' be considered for inclusion in the returned indices. If there are fewer 
#' TRUE values in \code{subset} than \code{n}, all of the unclassified data 
#' will be used as the subset.
#' @param n number of indices to return.
#'
#' @return An integer vector that contains the indices of the events that
#' should be included in the creation of the final classification model for
#' the parameter of interest (bead, debris, doublet, dead).
#'
#' @details
#' The indices that are returned by \code{modelData} should be used to
#' create a model that can be used to classify the observations with
#' regard to the parameter of interest (bead, doublet, debris, dead).
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' nn <- NN(x)
#' labels <- qcDataFrame(x)
#' beads <- beadID(x, labels = labels)
#' index <- modelData(labels, init = beads$init, n = 4000)
#'
#' @export
modelData <- function(labels, init, subset = NULL, n = 4000) {
    
    if (sum(init == 1) < (0.5 * n) | sum(init == -1) < (0.5 * n)) {
        warning(paste("Only", min(sum(init == 1), sum(init == -1)),
                      "observations in subset. All unclassified data used to select indices."))
    }
    
    # if (subset) {
    #     
    # }
    
    if (sum(init == 1) >= (0.5 * n) & sum(init == -1) >= (0.5 * n)) {
        poss.ind <- 1:(nrow(labels))
        poss.ind <- poss.ind[init != 0]
        poss.wt <- ifelse(init[poss.ind] == -1, (1000 / table(init[poss.ind]))[1], 
                          (1000 / table(init[poss.ind]))[2])
    } else {
        poss.ind <- which(labels$label == "cell")
        poss.wt <- ifelse(init[poss.ind] == -1, (1000 / table(init[poss.ind]))[1], 
                          (1000 / table(init[poss.ind]))[2])
    }
    
    sample(poss.ind, n, prob = poss.wt)
    
}
