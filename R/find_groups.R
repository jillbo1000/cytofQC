#' splits a numeric vector into groups.
#'
#' @param x A numeric vector to be split into groups.
#'
#' @return A vector that contains the groups assignment for each element
#' of \code{x}.
#'
#' @details
#' The function determines how many peaks are in the data and then
#' splits the data into groups that belong to each peak.
#'
#' @examples
#' fname <- "../data/FlowRepository_FR-FCM-Z29V_files/REP_1_deid.fcs"
#' x <- dataPrep(fname)
#' find_groups(x$Bead1)
#'
#' @export
find_groups <- function(x){

  d <- stats::density(x, adjust = 2)
  # initial peaks:
  #  (1) higher than it's 10 neighbors
  #  (2) at leat 1% of the maximum peak height
  peakTF <- sapply(1:length(d$y), function(i) {
    neighbors <- c((i - 5):(i - 1), (i + 1):(i + 5))
    neighbors <- neighbors[neighbors %in% 1:length(d$y)]

    return(all(d$y[i] > d$y[neighbors]) & d$y[i] > 0.01 * max(d$y))
  })

  peaks <- which(peakTF)
  # find valleys (lowpoints between peaks)
  if(sum(peakTF) > 1){
    valleys <- sapply(1:(sum(peakTF) - 1), function(i) {
      p1 <- which(peakTF)[i]
      p2 <- which(peakTF)[i + 1]
      btwn <- p1:p2
      return(btwn[which.min(d$y[btwn])])
    })
  }

  # peak filtering
  # remove peak if not >120% of both nearby valleys
  for (ii in 1:length(peaks)) {
    ht <- d$y[peaks[ii]]
    if(ii == 1){
      lv <- 1
    } else {
      lv <- valleys[ii - 1]
    }
    if (ii == length(peaks)) {
      rv <- length(d$x)
    } else {
      rv <- valleys[ii]
    }
    lv <- d$y[lv]
    rv <- d$y[rv]
    if(ht < 1.2 * lv | ht < 1.2 * rv){
      peakTF[peaks[ii]] <- FALSE
      print(paste('removing', ii))
    }
  }
  peaks <- which(peakTF)
  # re-find valleys
  if (sum(peakTF) > 1) {
    valleys <- sapply(1:(sum(peakTF) - 1), function(i){
      p1 <- which(peakTF)[i]
      p2 <- which(peakTF)[i + 1]
      btwn <- p1:p2
      return(btwn[which.min(d$y[btwn])])
    })
    out <- cut(x, breaks = c(-Inf, d$x[valleys], Inf))
  } else {
    out <- factor(rep(1, length(x)))
  }
  levels(out) <- 1:length(peaks)
  return(out)
}
