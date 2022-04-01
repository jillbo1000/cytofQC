#' Read in a dataset and prepare it for analysis
#'
#' @param file.name A path to an .fcs file that contains CyTOF data.
#' @param beads character vector that contains the names of all of the bead channels.
#' @param dna Character vector that contains the names of the DNA markers.
#' @param event_length Character vector of the event length variable.
#' @param viability Character vector of the permeability/viability markers.
#'
#' @return A \code{matrix} that contains the time, beads, DNA measures, permeability
#' (viability), event length, and the Guassian parameters of the CyTOF data.
#' Note that the time variable is included so that the cells can be matched up
#' exactly with the QC label data at any point in the process.
#'
#' @details
#' The function returns a \code{matrix} that contains the beads, DNA measures,
#' event length, permeability (viability) measure, and Gaussian parameters.
#' The event length and Gaussian parameters are transformed with a natural
#' log transformation. The data are not standardized at this point because
#' the raw data are needed for some of the data labeling steps.
#'
#' @examples
#' data("raw_data", package = "CATALYST")
#' tech <- dataPrep(raw_data, beads = 'Beads', viability = c('cisPt1','cisPt2'))
#' 
#' @export
dataPrep <- function(file.name,
                     time = "Time",
                     beads = c("Bead"),
                     dna = c("DNA1", "DNA2"),
                     event_length = "Event_length",
                     viability = "Live_Dead",
                     gaussian = c("Center", "Offset", "Width", "Residual")) {

  sce <- CATALYST::prepData(file.name)
  time_channel <- matrix(SingleCellExperiment::int_colData(sce)[, names(SingleCellExperiment::int_colData(sce)) %in% time],
                                                                nrow = ncol(SummarizedExperiment::assay(sce, "exprs")))
  colnames(time_channel) <- "Time"

  bead_channels <- matrix(t(SummarizedExperiment::assay(sce, "exprs")[rownames(sce) %in% beads, ]),
                          nrow = ncol(SummarizedExperiment::assay(sce, "exprs")))
  colnames(bead_channels) <- paste0("Bead", 1:ncol(bead_channels))

  dna_channels <- matrix(t(SummarizedExperiment::assay(sce, "exprs")[rownames(sce) %in% dna, ]),
                         nrow = ncol(SummarizedExperiment::assay(sce, "exprs")))
  if (ncol(dna_channels) > 1) {
    colnames(dna_channels) <- paste0("DNA", 1:ncol(dna_channels))
  } else {
    colnames(dna_channels) <- "DNA"
  }

  perm_channels <- matrix(t(SummarizedExperiment::assay(sce, "exprs")[rownames(sce) %in% viability, ]),
                          nrow = ncol(SummarizedExperiment::assay(sce, "exprs")))
  if (ncol(perm_channels) > 1) {
    colnames(perm_channels) <- paste0("Viability", 1:ncol(perm_channels))
  } else {
    colnames(perm_channels) <- "Viability"
  }

  gauss <- log1p(as.matrix(SingleCellExperiment::int_colData(sce)[, names(SingleCellExperiment::int_colData(sce)) %in% c(event_length, gaussian)]))

  # Make sure names of Gaussian variables are standardized
  colnames(gauss)[grep(event_length, colnames(gauss))] <- "Event_length"
  colnames(gauss)[grep("enter", colnames(gauss), ignore.case = TRUE)] <- "Center"
  colnames(gauss)[grep("off", colnames(gauss), ignore.case = TRUE)] <- "Offset"
  colnames(gauss)[grep("res", colnames(gauss), ignore.case = TRUE)] <- "Residual"
  colnames(gauss)[grep("wid", colnames(gauss), ignore.case = TRUE)] <- "Width"

  tech <- cbind(time_channel, bead_channels, dna_channels, perm_channels, gauss)
  tech
}
