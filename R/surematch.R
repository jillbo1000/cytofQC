#' Returns the match status for a parameter
#'
#' @param nn A matrix containing the indices of the nearest neighbors computed
#' with \code{\link{NN}}.
#' @param init A logical vector that contains the initial labeling for the
#' cells for the cell type of interest.
#' @param threshold Number of mismatches that classify an observation as a
#' not a sure match. The default is set to 0, which means that an observation
#' is classified as a match only if it matches all of its neighbors. 
#'
#' @return A \code{matrix} that contains the match status for each
#' observation and its neighbors. A value of TRUE in the matrix indicates
#' that the neighbor is mismatched with the observation for the row. The
#' last column is the match status for the observation. If the sum of
#' TRUE values in the row is less than or equal to the threshold, the
#' observation is classified as a match, indicated with a value of TRUE.
#' These TRUE values are used to select a set of cells that are candidates
#' for training the classification model.
#'
#' @details
#' This function uses the the nearest neighbors and the initial labels of the
#' unclassified cells to identify a set of cells that are different than their
#' neighbors in regard to its labeling. This is done by giving each index in
#' NN matrix a TRUE value if \code{init} is true for that observation and a
#' FALSE if \code{init} is FALSE. If the observation itself has a value of
#' TRUE in the init vector, the values in its row are flipped. This
#' creates a matrix that has a value of TRUE if the observation and its
#' neighbor are mismatched. NA values in the NN matrix are all set to FALSE.
#' An observation is classified as a mismatch if \code{m} or more of its
#' neighbors are a mismatch (denoted by TRUE).
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' nn <- NN(x)
#' labels <- qcDataFrame(x)
#' beads <- initialBead(x, labels = labels)
#' mm <- surematch(nn, init = beads$init, threshold = 1)
#'
#'@export
surematch <- function(nn, init, threshold = 0) {

  Time <- nn[, 1]
  nn <- nn[, -1]

  mm <- matrix(init[nn], ncol = ncol(nn))
  mm[init, ] <- !mm[init, ]
  mm[is.na(mm)] <- FALSE
  colnames(mm) <- paste0("nn", 1:ncol(mm))

  status <- rowSums(mm) <= threshold

  mism <- data.frame(Time, mm, match = status)
  mism

}
