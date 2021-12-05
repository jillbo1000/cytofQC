#' Preliminary viability classification.
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}. Note that the
#' data should not be standardized. All of the bead channels should be
#' non-negative.
#' @param labels A \code{data.frame} created with \code{\link{qcDataFrame}}.
#' @param dna If TRUE, the DNA will be used to determine viability.
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' this analysis.
#'
#' @return A \code{data.frame} that contains the viablity designation
#' for the observation. The first column is the \code{Time} variable
#' that helps ensure that all of the observations can be correctly matched
#' during all stages labeling. The last column is \code{init} which is
#' TRUE if the observation is determined to be dead and FALSE if the
#' observation is not dead during this initial classification.
#'
#' @details
#' The beads are typically the first cell classification that is done. The
#' different cell types are labeled iteratively so the \code{labels}
#' data.frame should contain all of the labels and probabilities computed
#' up to this point. Thus, if the dead cells were not the first cell type that
#' was classified, the \code{labels} data.frame must have the labels and
#' probabilities computed prior to classifying the dead cells. Also note that
#' this matrix will not be changed during this state. This function
#' creates the most preliminary bead labels and that information is needed
#' to obtain the final bead labels. The \code{labels} data.frame is needed
#' to determine which cells have and have not been labeled prior to this
#' step.
#'
#' The viability measures should classify into two fairly clear groups
#' where one is permeable cells and the other is non-permeability. DNA is
#' also often higher for cells that are permeable. The primary measure
#' for determining permeability is sum of the viability and the DNA measures.
#' The function \code{find_groups} is used to determine the groups. The members of the group with the
#' largest mean are classified as beads and the rest are classified as
#' not beads. The observations that are labeled 'GDPzero' are included in
#' the output data, but the are removed from the data for bead designation.
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' labels <- qcDataFrame(x)
#' dead <- initialDead(x, labels = labels)
#'
#' @export
initialDead <- function(x, labels, dna = TRUE, standardize = TRUE) {

  Time <- x[, 1]
  x <- x[, -1]

  if (standardize) {
    x <- scale(x)
  }

  unclassified.ind <- which(labels$label == "cell")
  cell <- x[unclassified.ind, ]

  viability <- x[, grep("Viability", colnames(x))]

  if (ncol(data.frame(viability)) > 1) {
    viability <- rowSums(viability, na.rm = TRUE)
  }

  if (dna) {
    dnaData <- x[, grep("DNA", colnames(x))]
    if (ncol(dnaData) > 1) {
      n <- ncol(dnaData)
      dnaData <- rowSums(dnaData, na.rm = TRUE)
      dnaData <- dnaData / n
    }
    deadScore <- viability + 0.5 * dnaData
  } else {
    deadScore <- viability
  }

  live.ind <- which(viability == min(viability))
  live.un <- live.ind[live.ind %in% unclassified.ind]

  g <- find_groups(deadScore[live.un])
  mns <- by(deadScore[live.un], g, mean)
  deadclus <- which.max(mns)
  liveclus <- which.min(mns)

  gAll <- rep(NA, nrow(x))
  gAll[unclassified.ind] <- deadclus
  gAll[live.un] <- liveclus

  init <- rep(FALSE, nrow(x))

  if (length(unique(g)) > 1) {
    init[unclassified.ind] <- gAll[unclassified.ind] == deadclus
  }

  data.frame(Time, deadScore = deadScore, group = gAll, init = init)

}
