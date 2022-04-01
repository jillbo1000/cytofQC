#' Creates a \code{matrix} that contains the nearest neighbors indices
#'
#' @param x A \code{matrix} created with \code{\link{dataPrep}}.
#' @param n Number of neighbors to include in the matrix.
#' @param type "knn" will do k-nearest neighbors and "mnn" will do mutual
#' nearest neighbors.
#' @param standardize Indicates if the data should be standardized. Because
#' the data are on different scales, it should be standardized for
#' any analysis.
#'
#' @return A \code{matrix} that contains the indices of the nearest neighbors
#' for all of the observations in \code{x}. If "mnn" is selected, there will
#' be many NAs in the matrix which indicate neighbors for an observation that
#' are not mutual.
#'
#' @details
#' This matrix is used prior to computing the classification models for each
#' cell type. It is used to identify mismatched cells and only needs to be
#' computed once.
#' NN is done using all of the QC parameters. The time variable is removed and
#' the remaining variables are standardized. NN is done with the standardized
#' data. The time variable is included with the matrix that is returned so
#' that the order of the data can be preserved for all data objects and all
#' steps of the process.
#'
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#'
#' # See names of variables
#' require(CATALYST)
#' require(SingleCellExperiment)
#' sce <- prepData(fname)
#' rownames(sce)
#' rowData(sce)
#' names(int_colData(sce))
#'
#' # Read in data
#' x <- dataPrep(fname)
#' knn <- NN(x, n = 10)
#' mnn <- NN(x, )
#'
#' @export
NN <- function(x, n = 20, type = "mnn", standardize = TRUE) {

  Time <- x[, 1]
  x <- x[, -1]

  if(standardize) {
    x <- scale(x)
  }

  k.dex <- BiocNeighbors::buildKmknn(x)
  knn <- BiocNeighbors::findKNN(x, n, BNINDEX = k.dex)

  if (type == "mnn") {
    mnn <- knn$index
    for (i in 1:nrow(mnn)) {
      index <- knn$index[i, ]
      keep <- rowMeans(knn$index[index, ] == i) > 0
      mnn[i, !keep] <- NA
    }
    cbind(Time, mnn)
  } else {
    cbind(Time, knn$index)
  }

}
