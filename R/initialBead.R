#' Preliminary bead classification.
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}. Note that the
#' data should not be standardized. All of the bead channels should be
#' non-negative.
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#'
#' @return A \code{data.frame} that contains the bead designation for each
#' bead channel and the overall bead designation for the observation. The
#' first column is the \code{Time} variable that helps ensure that all
#' of the observations can be correctly matched during all stages labeling.
#' The last column is \code{init} which is TRUE if the observation is
#' determined to be a bead and FALSE if the observation is not a bead during
#' this initial classification.
#'
#' @details
#' The beads are typically the first cell classification that is done. The
#' different cell types are labeled iteratively so the \code{labels}
#' data.frame should contain all of the labels and probabilities computed
#' up to this point. Thus, if the beads were not the first cell type that
#' was classified, the \code{labels} data.frame must have the labels and
#' probabilities computed prior to classifying the beads. Also note that
#' this matrix will not be changed during this state. This function
#' creates the most preliminary bead labels and that information is needed
#' to obtain the final bead labels. The \code{labels} data.frame is needed
#' to determine which cells have and have not been labeled prior to this
#' step.
#'
#' Each bead channel should classify into two fairly clear groups where one
#' is the beads and the other is non-beads. The function \code{find_groups}
#' is used to determine the groups. The members of the group with the
#' largest mean are classified as beads and the rest are classified as
#' not beads. The observations that are labeled 'GDPzero' are included in
#' the output data, but the are removed from the data for bead designation.
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' labels <- qcDataFrame(x)
#' beads <- initialBead(x, labels = labels)
#'
#' @export
initialBead <- function(x, labels) {

  unclassified.ind <- which(labels$label == "cell")
  cell <- x[unclassified.ind, ]
  bead_channels <- grep("Bead", colnames(cell))

  if (min(bead_channels) < 0) {
    stop("Bead data should all be non-negative")
  }
  
  if (length(bead_channels) > 1) {
    b <- rowSums(x[, bead_channels], na.rm = TRUE)
  } else if (length(bead_channels) == 1){
    b <- x[, bead_channels]
  } else {
    stop("No bead channels in data")
  }
  
  g <- initialGuess(b[unclassified.ind])
  init <- rep(FALSE, nrow(x))
  init[unclassified.ind] <- (g$label == max(g$label))
  
  data.frame(Time = x[, "Time"], beadScore = b, init = init)
}
