% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/initialBead.R
\name{initialBead}
\alias{initialBead}
\title{Preliminary bead classification.}
\usage{
initialBead(x, labels, threshold = 0.5)
}
\arguments{
\item{x}{A \code{matrix} created with \code{\link{dataPrep}}. Note that the
data should not be standardized. All of the bead channels should be
non-negative.}

\item{labels}{A \code{data.frame} created with \code{\link{qcDataFrame}}.}

\item{threshold}{Number in [0, 1] that indicates the number of
bead channels an observation needs to be flagged as a bead to be
considered a bead.}
}
\value{
A \code{data.frame} that contains the bead designation for each
bead channel and the overall bead designation for the observation. The
first column is the \code{Time} variable that helps ensure that all
of the observations can be correctly matched during all stages labeling.
The last column is \code{init} which is TRUE if the observation is
determined to be a bead and FALSE if the observation is not a bead during
this initial classification.
}
\description{
Preliminary bead classification.
}
\details{
The beads are typically the first cell classification that is done. The
different cell types are labeled iteratively so the \code{labels}
data.frame should contain all of the labels and probabilities computed
up to this point. Thus, if the beads were not the first cell type that
was classified, the \code{labels} data.frame must have the labels and
probabilities computed prior to classifying the beads. Also note that
this matrix will not be changed during this state. This function
creates the most preliminary bead labels and that information is needed
to obtain the final bead labels. The \code{labels} data.frame is needed
to determine which cells have and have not been labeled prior to this
step.

Each bead channel should classify into two fairly clear groups where one
is the beads and the other is non-beads. The function \code{find_groups}
is used to determine the groups. The members of the group with the
largest mean are classified as beads and the rest are classified as
not beads. The observations that are labeled 'GDPzero' are included in
the output data, but the are removed from the data for bead designation.
}
\examples{
fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
x <- dataPrep(fname)
labels <- qcDataFrame(x)
beads <- initialBead(x, labels = labels)

}
