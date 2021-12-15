#' Preliminary doublet classification.
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}.
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param score Integer that identifies the doublet score that should be used.
#' See \code{Details} for a description of the different scores that can be
#' used.
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' this analysis.
#'
#' @return A \code{data.frame} that contains the doublet score and the
#' overall doublet designation for the observation. The variables included
#' the \code{data.frame} are:
#' \item{Time}{The time stamp from the original data. It is needed to ensure
#' that the labels can be matched with the original cell data.}
#' \item{doubletScore}{The computed doublet score.}
#' \item{doubletGroup}{The group that the score was assigned to. This is
#' different than the label because the doublets typically cluster into
#' three groups.}
#' \item{init}{TRUE if the observation is determined to be a doublet and
#' FALSE if the observation is not a doublet during this initial
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
#' up to this point. Thus, if the doublets were not the first cell type that
#' was classified, the \code{labels} data.frame must have the labels and
#' probabilities computed prior to classifying the doublets. Also note that
#' this matrix will not be changed during this state. This function
#' creates the initial doublet labels and that information is needed
#' to obtain the final doublet labels. The \code{labels} data.frame is needed
#' to determine which cells have and have not been labeled prior to this
#' step.
#'
#' Several options are available for computing the doublet score. The
#' following list shows the doublet score calculations. Each one can
#' be selected by its number on the following list:
#' \enumerate{
#'   \item{DNA1 + DNA2 + Residual + Event_length - Offset - 0.5(Width)}
#'   \item{DNA1 + DNA2 + Residual + Event_length - Offset - 0.5(Width) + abs(Center)}
#'   \item{DNA1 + DNA2 + 0.5(Residual + Event_length - Offset - 0.5(Width) + abs(Center))}
#' }
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' labels <- qcDataFrame(x)
#' doublets <- initialDebris(x, labels = labels, score = 3)
#'
#' @export
initialDoublet <- function(x, labels, score = 3, standardize = TRUE) {

  Time <- subset(x, select = c(get("Time")))
  x <- subset(x, select = -c(get("Time")))
  
  if (standardize) {
    x <- scale(x)
  }

  unclassified.ind <- which(labels$label == "cell")
  cell <- data.frame(x[unclassified.ind, ])

  if (score == 1) {
    doubletScore <- x[, "DNA1"] + x[, "DNA2"] + x[, "Residual"] +
      x[, "Event_length"] - x[, "Offset"] - 0.5 * x[, "Width"]
  } else if (score == 2) {
    doubletScore <- x[, "DNA1"] + x[, "DNA2"] + x[, "Residual"] +
      x[, "Event_length"] - x[, "Offset"] - 0.5 * x[, "Width"] +
      abs(x[, "Center"])
  } else if (score == 3) {
    doubletScore <- x[, "DNA1"] + x[, "DNA2"] +
      0.5 * (x[, "Residual"] + x[, "Event_length"] - x[, "Offset"]
             - 0.5 * x[, "Width"] + abs(x[, "Center"]))
  } else {
    stop("Invalid score selection")
  }

  g <- find_groups(doubletScore[unclassified.ind])
  mns <- by(doubletScore[unclassified.ind], g, mean)
  doubletclus <- which.max(mns)
  init <- rep(FALSE, nrow(x))
  init[unclassified.ind] <- g == doubletclus

  groups <- rep(NA, nrow(x))
  groups[unclassified.ind] <- g

  data.frame(Time, doubletScore = doubletScore, doubletGroup = groups, init = init)

}
