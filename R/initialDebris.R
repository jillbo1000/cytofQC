#' Preliminary debris classification.
#'
#' @import mclust
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}.
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param score Integer that identifies the debris score that should be used.
#' See \code{Details} for a description of the different scores that can be
#' used.
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' this analysis.
#'
#' @return A \code{data.frame} that contains the debris score and the
#' overall debris designation for the observation. The variables included
#' the \code{data.frame} are:
#' \item{Time}{The time stamp from the original data. It is needed to ensure
#' that the labels can be matched with the original cell data.}
#' \item{debrisScore}{The computed debris score.}
#' \item{debrisGroup}{The group that the score was assigned to. This is
#' may be different than the label because the debris often cluster into
#' two groups.}
#' \item{init}{TRUE if the observation is determined to be debris and
#' FALSE if the observation is not debris during this initial
#' classification.}
#'
#' @details
#' The beads are typically the first cell classification that is done
#' because their identification is straightforward. Doublets and
#' debris are more difficult to identify so they should be classified
#' after the beads.
#'
#' The different cell types are labeled iteratively so the \code{labels}
#' data.frame should contain all of the labels and probabilities computed
#' up to this point. Thus, if the debris were not the first cell type that
#' was classified, the \code{labels} data.frame must have the labels and
#' probabilities computed prior to classifying the debris. Also note that
#' this matrix will not be changed during this state. This function
#' creates the initial debris labels and that information is needed
#' to obtain the final debris labels. The \code{labels} data.frame is needed
#' to determine which cells have and have not been labeled prior to this
#' step.
#'
#' Several options are available for computing the debris score. The
#' following list shows the debris score calculations. Each one can
#' be selected by its number on the following list:
#' \enumerate{
#'   \item{1 - (DNA1 + DNA2 + Event_length - Center - Width + Offset)}
#'   \item{Residual + Offset - 2(DNA1) - 2(DNA2) - Event_length - Center - 0.5(Width)}
#'   \item{1 - (DNA1 + DNA2 + Event_length)}
#' }
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' labels <- qcDataFrame(x)
#' debris <- initialDebris(x, labels = labels, score = 3)
#'
#' @export
initialDebris <- function(x, labels, score = 3, standardize = TRUE) {

  if (standardize) {
    xs <- scale(x[, -1])
  } else {
    xs <- x
  }
  
  unclassified.ind <- which(labels$label == "cell")
  cell <- data.frame(xs[unclassified.ind, ])
  
  if (score == 1) {
    debrisScore <- 1 - (x[, "DNA1"] + x[, "DNA2"] + x[, "Event_length"] -
                          x[, "Center"] - x[, "Width"] + x[, "Offset"])
  } else if (score == 2) {
    debrisScore <- x[, "Residual"] + x[, "Offset"] - 2.0 * x[, "DNA1"] -
      2.0 * x[, "DNA2"] - x[, "Event_length"] - x[, "Center"] - 0.5 * x[, "Width"]
  } else if (score == 3) {
    debrisScore <- 1 - (x[, "DNA1"] + x[, "DNA2"] + x[, "Event_length"])
  } else {
    stop("Invalid score selection")
  }
  
  g <- initialGuess(debrisScore[unclassified.ind], middleGroup = 1)
  init <- rep(0, nrow(x))
  init[unclassified.ind] <- (g$label != min(g$label))
  
  data.frame(Time = x[, 1], debrisScore = debrisScore, init = init)

}
