#' @importFrom graphics hist abline curve layout
NULL
#' Plot preliminary classification from initialGuess
#'
#' @param x The score (ie. debris score, doublet score, etc.) to be used for
#'   predicting each event's label (eg. "doublet" vs. "cell").
#' @param IG If NULL, the function \code{initialGuess} is used to fit a mixture
#'   of normal distributions. Otherwise, a numeric vector can be passed to the 
#'   function that contains the fitted values.
#' @param fit If left blank or NULL, the best fit as determined by BIC will be
#'   plotted. Otherwise, a numeric value of 1, 2, or 3 can be be selected to 
#'   plot single normal fit, mixture of two normals, or the mixture of three
#'   normals as fit by \code{initialGuess}. 
#' @param type Type of graph to be plotted. If 'truncated' is selected, only 
#'   half of the first normal distribution will be plotted. If 'both' is 
#'   selected, both the truncated and full plot will be plotted. 
#' 
#' @return A histogram that shows the score with the mixture of normal 
#' distributions overlayed. 
#'   
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- initialDoublet(sce)
#' plotInitialGuess(scores(sce, "doublet"), type = "both")
#' plotInitialGuess(scores(sce, "doublet"), type = "truncated")
#' 
#' @import graphics
#' @export

plotInitialGuess <- function(x, IG = NULL, fit = NULL, 
                             type = c("both", "full", "truncated")){
    
    if (!methods::is(x, "numeric")) {
        stop("x must be a numeric vector")
    }
    
    type <- match.arg(tolower(type), 
                      choices = c("both", "full", "truncated"))
    
    
    if(is.null(IG)){
        IG <- initialGuess(x)
    }
    
    d <- density(x[which(x > min(x))]) 
    cut <- d$x[which.max(d$y)]
    xx <- x[x >= cut]
    xx <- xx - cut
    p0 <- mean(x >= cut)
    
    if(is.null(fit)){
        fit <- which.min(c(IG$fit1$bic, IG$fit2$bic, IG$fit3$bic))
    }
    
    if (type == "both") {
        layout(matrix(c(1, 2), nrow=1))
    }
    
    if (type %in% c("both", "full")) {
        hist(x, breaks=100, probability = TRUE, main = 'Full Score')
        abline(v = cut)
        if(fit == 1){
            curve(p0*2*dnorm(x-cut,sd=IG$fit1$pars[1]), from = cut, 
                  to = max(x)+1, col = 4, lwd=2, add = TRUE)
        }
        if(fit == 2){
            p1 <- IG$fit2$pars[1]
            s1 <- IG$fit2$pars[2]
            m2 <- IG$fit2$pars[3]
            s2 <- IG$fit2$pars[4]
            curve(p0*2*p1*dnorm(x-cut,sd=s1), from = cut, to = max(x)+1, 
                  col = 4, lwd=2, add = TRUE)
            curve(p0*(1-p1)*dnorm(x-cut,mean=m2,sd=s2), from = cut, 
                  to = max(x)+1, col = 2, lwd=2, add = TRUE)
        }
        if(fit == 3){
            p1 <- IG$fit3$pars[1]; p2 <- IG$fit3$pars[2]; p3 <- IG$fit3$pars[3]
            m2 <- IG$fit3$pars[4]; m3 <- IG$fit3$pars[5]
            s1 <- IG$fit3$pars[6]; s2 <- IG$fit3$pars[7]; s3 <- IG$fit3$pars[8]
            curve(p0*2*p1*dnorm(x-cut,sd=s1), from = cut, to = max(x)+1, col = 4, 
                  lwd=2, add = TRUE)
            curve(p0*p2*dnorm(x-cut,mean=m2,sd=s2), from = cut, to = max(x)+1, 
                  col = 1, lwd=2, add = TRUE)
            curve(p0*p3*dnorm(x-cut,mean=m3,sd=s3), from = cut, to = max(x)+1, 
                  col = 2, lwd=2, add = TRUE)
        }
    }
    
    if (type %in% c("both", "truncated")) {
        hist(xx, breaks=100, probability = TRUE, main = 'Truncated Score')
        if(fit == 1){
            curve(2*dnorm(x,sd=IG$fit1$pars[1]), from = 0, to = max(x)+1, col = 4, 
                  lwd=2, add = TRUE)
        }
        if(fit == 2){
            p1 <- IG$fit2$pars[1]
            s1 <- IG$fit2$pars[2]
            m2 <- IG$fit2$pars[3]
            s2 <- IG$fit2$pars[4]
            curve(2*p1*dnorm(x,sd=s1), from = 0, to = max(xx)+1, col = 4, 
                  lwd=2, add = TRUE)
            curve((1-p1)*dnorm(x,mean=m2,sd=s2), from = 0, to = max(xx)+1, 
                  col = 2, lwd=2, add = TRUE)
        }
        if(fit == 3){
            p1 <- IG$fit3$pars[1]; p2 <- IG$fit3$pars[2]; p3 <- IG$fit3$pars[3]
            m2 <- IG$fit3$pars[4]; m3 <- IG$fit3$pars[5]
            s1 <- IG$fit3$pars[6]; s2 <- IG$fit3$pars[7]; s3 <- IG$fit3$pars[8]
            curve(2*p1*dnorm(x,sd=s1), from = 0, to = max(xx)+1, col = 4, lwd=2, 
                  add = TRUE)
            curve(p2*dnorm(x,mean=m2,sd=s2), from = 0, to = max(xx)+1, col = 1, 
                  lwd=2, add = TRUE)
            curve(p3*dnorm(x,mean=m3,sd=s3), from = 0, to = max(xx)+1, col = 2, 
                  lwd=2, add = TRUE)
        }
        
    }
    layout(1)
}

