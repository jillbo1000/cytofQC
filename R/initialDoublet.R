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
#'   \item{0.3 * (DNA1 + DNA2 + Event_length) + Residual + Center + (max(Offset))}
#' }
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' tech <- dataPrep(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' labels <- qcDataFrame(tech)
#' initialDebris(tech, labels = labels, score = 3)
#'
#' @export
initialDoublet <- function(x, labels, score = 3, standardize = TRUE) {

  if (standardize) {
    xs <- scale(x[, -1])
  } else {
    xs <- x
  }

  unclassified.ind <- which(labels$label == "cell")
  cell <- data.frame(xs[unclassified.ind, ])

  if (score == 1) {
    doubletScore <- xs[, "DNA1"] + xs[, "DNA2"] + xs[, "Residual"] +
      xs[, "Event_length"] - xs[, "Offset"] - 0.5 * xs[, "Width"]
  } else if (score == 2) {
    doubletScore <- xs[, "DNA1"] + xs[, "DNA2"] + xs[, "Residual"] +
      xs[, "Event_length"] - xs[, "Offset"] - 0.5 * xs[, "Width"] +
      abs(xs[, "Center"])
  } else if (score == 3) {
    doubletScore <- 0.3 * (xs[, "DNA1"] + xs[, "DNA2"] +
      xs[, "Event_length"]) + xs[, "Residual"] + xs[, "Center"] +
      (max(xs[, "Offset"]) - xs[, "Offset"])
  } else {
    stop("Invalid score selection")
  }

  g <- initialGuess(doubletScore[unclassified.ind], middleGroup = 1)
  init <- rep(0, nrow(x))
  init[unclassified.ind] <- g$label
  
  data.frame(Time = x[, 1], doubletScore = doubletScore, init = init)
}
