#' Returns indices for data to be used to create the final classification model.
#'
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param subset A logical vector that indicates if an observation should
#' be considered for inclusion in the returned indices. This vector may
#' be computed with \code{\link{mismatch}}. If there are fewer TRUE values
#' in \code{subset} than \code{n}, all of the unclassified data will
#' be used as the subset.
#' @param init A logical vector that contains the initial labeling for the
#' observations for the parameter of interest.
#' @param n number of indices to return.
#'
#' @return An integer vector that contains the indices of the datapoints that
#' should be included in the creation of the final classification model for
#' the parameter of interest (bead, doublet, debris, dead).
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
#' mm <- mismatch(nn, init = beads$init, threshold = 1)
#' index <- modelData(labels, subset = mm$Mismatch, init = beads$init, n = 4000)
#'
#' @export
modelData <- function(labels, subset, init, n = 4000) {

  if (sum(subset) < n) {
    warning(paste("Only", sum(subset),
                  "observations in subset. All unclassified data used to select indices."))
  }

  if (sum(subset) >= n) {
    poss.ind <- 1:(nrow(labels))
    poss.ind <- poss.ind[subset]
    poss.wt <- (1000 / table(init[poss.ind]))[1 + init[poss.ind]]
  } else {
    poss.ind <- which(labels$label == "cell")
    poss.wt <- (1000 / table(init[poss.ind]))[1 + init[poss.ind]]
  }

  sample(poss.ind, n, prob = poss.wt)

}
