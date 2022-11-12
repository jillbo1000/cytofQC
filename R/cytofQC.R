#' cytofQC: Event labeling for quality control of CyTOF data
#'
#' Labels observations in a CyTOF dataset as a cell, gdpZero (zero for at least 
#' one Gaussian parameter), bead, debris, doublet, or dead cell.  
#' 
#'  
#' Data from an fcs file are read directly into a 
#' \code{\link{SingleCellExperiment}} using the \code{\link{readCytof}} 
#' function. 
#' 
#' The data can be labeled with a single function, \code{\link{labelQC}}, which
#' can be customized. Labeling can also be done using a set of other functions
#' that first select a set of events that clearly look like the event type 
#' being modeled and then use those events to train a statistical learning 
#' model that can identify the event type. These functions are discussed and
#' demonstrated in the vignette. 
#' 
#' A plotting function called \code{\link{cytofHist}} is included that makes 
#' assessing the characteristic of the data and labeling easy. 
#' 
#' The package also includes a function called \code{\link{cytofQCreport}} 
#' that generates a report of the labeling and can generate a umap created with 
#' the QC variables and colored by event label. 
#' 
#' \tabular{ll}{ Package: \tab cytofQC \cr Type: \tab Package \cr Version:
#' \tab 0.99.3 \cr Date: \tab 2022-11-11 \cr License: \tab Artistic-2.0\cr }
#' 
#' @name cytofQC-package
#' @aliases cytofQC cytofQC-package
#' @docType package
#' @author 
#' Maintainer: Jill Lundell <jflundell@gmail.com>
#' 
#' Authors: J. Lundell, K. Street
#' 
NULL