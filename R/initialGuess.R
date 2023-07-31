# produce initial guess labels based on mixture model 
# (truncated normal + normals)

# fit a half normal
#' @import stats
.mixFit1 <- function(x){
    s0 <- sqrt(sum(x^2) / (length(x)-1))
    psingle <- dnorm(x, sd = s0, log = TRUE) + log(2)
    ll <- sum(psingle)
    return(list(pars = c(s0=s0), ll = ll, bic = -2*ll+1*log(length(x)),
                dens = matrix(exp(psingle))))
}

# fit a mixture of a half normal and a normal
#' @importFrom matrixStats weightedSd
.mixFit2 <- function(x, thresh = .01){
    # hidden factor, Z in [0,1], is prob observation X originated from 
    # left distn.
    # p1 := mean(Z)
    Z <- as.numeric(x <= max(x)*1/2)
    if(sum(1-Z) < 3){
        Z <- as.numeric(x < sort(x, decreasing = TRUE)[3])
    }
    
    Z.old <- rep(1, length(x))
    p1 <- mean(Z)
    m2 <- weighted.mean(x, w = 1-Z) 
    s1 <- sqrt(sum(Z*(x^2)) / (sum(Z)-1))
    s2 <- weightedSd(x, 1-Z)
    
    working <- TRUE
    trace <- NULL
    while(working){
        pleft <- dnorm(x, mean = 0, sd = s1) * p1 * 2
        #pleft <- dt(x/s1, df = max(c(1,length(x)*p1-1))) * p1 * 2
        pright <- dnorm(x, mean = m2, sd = s2) * (1-p1)
        S <- pleft + pright
        Z <- pleft / S
        Z[S==0] <- .5
        
        if(max(abs(Z - Z.old)) < thresh){
            working <- FALSE
        }
        
        trace <- c(trace, sum(abs(Z-Z.old)))
        Z.old <- Z
        p1 <- mean(Z)
        
        s1.old <- s1
        s2.old <- s2
        
        m2 <- weighted.mean(x, w = 1-Z) 
        s1 <- sqrt(sum(Z*(x^2)) / (sum(Z)-1))
        s2 <- weightedSd(x, 1-Z)
        
        if(is.na(s1) | s1 == 0) s1 <- s1.old
        if(is.na(s2) | s2 == 0) s2 <- s2.old
    }
    pleft <- dnorm(x, mean = 0, sd = s1) * p1 * 2
    pright <- dnorm(x, mean = m2, sd = s2) * (1-p1)
    S <- pleft + pright
    S[S==0] <- .Machine$double.eps
    ll <- sum(log(S))
    
    return(list(pars = c(p1=p1,s1=s1,m2=m2,s2=s2), 
                ll = ll, 
                bic = -2*ll+4*log(length(x)),
                dens = cbind(pleft, pright)))
}

# fit a mixture of a half normal and two normals
.mixFit3 <- function(x, thresh = .01){
    # 3-component mixture model
    # hidden factor, Z1 in [0,1], prob observation X originated from left distn.
    # Z2 in [0,1] is middle distn prob
    # Z3 = 1-(Z1+Z2)
    # p1 := mean(Z1)
    # p2 := mean(Z2)
    Z <- Z.old <- matrix(0, nrow = length(x), ncol = 3)
    Z[,1] <- as.numeric(x <= max(x)/3)
    Z[,3] <- as.numeric(x > max(x)*2/3)
    if(sum(Z[,3]) < 3){
        Z[,3] <- as.numeric(x >= sort(x, decreasing = TRUE)[3])
    }
    Z[,2] <- 1 - Z[,1] - Z[,3]
    
    p <- colMeans(Z)
    
    m2 <- weighted.mean(x, w = Z[,2])
    m3 <- weighted.mean(x, w = Z[,3])
    
    s1 <- sqrt(sum(Z[,1]*(x^2)) / (sum(Z[,1])-1))
    s2 <- weightedSd(x, Z[,2])
    s3 <- weightedSd(x, Z[,3])
    
    working <- TRUE
    trace <- NULL
    while(working){
        pleft <- dnorm(x, mean = 0, sd = s1) * p[1] * 2
        pmid <- dnorm(x, mean = m2, sd = s2) * p[2]
        pright <- dnorm(x, mean = m3, sd = s3) * p[3]
        S <- pleft + pmid + pright
        
        Z <- cbind(pleft, pmid, pright) / S
        Z[S==0, ] <- 1/3
        
        if(max(abs(Z - Z.old)) < thresh){
            working <- FALSE
        }
        
        trace <- c(trace, sum(abs(Z-Z.old)))
        Z.old <- Z
        
        p <- unname(colMeans(Z))
        
        s1.old <- s1
        s2.old <- s2
        s3.old <- s3
        # m2.old <- m2
        # m3.old <- m3
        
        m2 <- weighted.mean(x, w = Z[,2])
        m3 <- weighted.mean(x, w = Z[,3])
        
        suppressWarnings(s1 <- sqrt(sum(Z[,1]*(x^2)) / (sum(Z[,1])-1)))
        suppressWarnings(s2 <- weightedSd(x, Z[,2]))
        suppressWarnings(s3 <- weightedSd(x, Z[,3], na.rm = TRUE))
        
        if(is.na(s1) | s1 == 0) s1 <- s1.old
        if(is.na(s2) | s2 == 0) s2 <- s2.old
        if(is.na(s3) | s3 == 0) s3 <- s3.old
    }
    # ensure group 3 is highest
    if(m2 > m3){
        m4 <- m3
        s4 <- s3
        p4 <- p[3]
        
        m3 <- m2
        s3 <- s2
        p[3] <- p[2]
        
        m2 <- m4
        s2 <- s4
        p[2] <- p4
    }
    pleft <- dnorm(x, mean = 0, sd = s1) * p[1] * 2
    pmid <- dnorm(x, mean = m2, sd = s2) * p[2]
    pright <- dnorm(x, mean = m3, sd = s3) * p[3]
    S <- pleft + pmid + pright
    S[S==0] <- .Machine$double.eps
    ll <- sum(log(S))
    
    return(list(pars = c(p1=p[1],p2=p[2],p3=p[3], 
                         m2=m2,m3=m3, s1=s1,s2=s2,s3=s3), 
                ll = ll, bic = -2*ll+7*log(length(x)),
                dens = cbind(pleft,pmid,pright)))
}


#' General preliminary classification.
#'
#' @param x The score (ie. debris score, doublet score, etc.) to be used for
#'   predicting each event's label (eg. "doublet" vs. "cell").
#' @param middleGroup numeric. When the optimal model (according to BIC) is the
#'   3-component mixture model, this argument determines how to assign the
#'   middle group. Possible values are \code{-1} for "cell", \code{0} (default)
#'   for "indeterminate", and \code{1} for the event type of interest (eg.
#'   "doublet").
#' @param bead logical. Should be TRUE when classifying beads. The bead score 
#'   sometimes has a larger peak for for the group that contains the beads. 
#'   This results in a misclassification of the beads that can make the 
#'   algorithm unable to find the beads. An adjustment is made to ensure this
#'   does not happen, but it should only be applied for bead classification.
#'   If the correction is used for other event types, they may be 
#'   misclassified.  
#' 
#' @return A list with the following elements: \itemize{ \item{\code{label}} {A
#'   vector of the same length as \code{x} providing the labels (\code{-1} for
#'   cells, \code{1} for non-cells, \code{0} for uncertain).} \item{\code{fit1}}
#'   {Summary of the 1-component (half Normal) model fit.} \item{\code{fit2}}
#'   {Summary of the 2-component (half Normal + Normal) model fit.}
#'   \item{\code{fit1}} {Summary of the 3-component (half Normal + 2 Normals)
#'   model fit.}}
#'   
#' @examples
#' data("raw_data", package = "CATALYST")
#' sce <- readCytof(raw_data, beads = "Beads", viability = c("cisPt1", "cisPt2"))
#' sce <- initialBead(sce)
#' fit <- initialGuess(scores(sce, "bead"), bead = TRUE)
#' sce <- initialDoublet(sce)
#' fit <- initialGuess(scores(sce, "doublet"))
#' 
#' @export
initialGuess <- function(x, middleGroup = c(0, -1, 1), bead = FALSE){
    
    if (!methods::is(x, "numeric")) {
        stop("x must be a numeric vector")
    }
    
    middleGroup <- as.numeric(match.arg(as.character(middleGroup), 
                                        c("0", "-1", "1")))

    if (bead) {
        mm <- mixtools::normalmixEM(x[x > min(x, na.rm = TRUE)])
        cut <- min(mm$mu)
    } else {
        d <- density(x[which(x > min(x))]) 
        cut <- d$x[which.max(d$y)]
    }

    xx <- x[x >= cut]
    xx <- xx - cut
    
    fit1 <- .mixFit1(xx)
    fit2 <- .mixFit2(xx)
    fit3 <- .mixFit3(xx)
    
    which.fit <- which.min(c(fit1$bic, fit2$bic, fit3$bic))
    
    if(which.fit == 1){
        # pick one group
        lab <- rep(-1, length(x))
    }else if (which.fit == 2){
        # pick two groups
        lab <- rep(-1, length(x))
        labxx <- rep(0, length(xx))
        labxx[fit2$dens[,1] > 99*fit2$dens[,2]] <- -1
        labxx[fit2$dens[,2] > 99*fit2$dens[,1]] <- 1
        lab[x >= cut] <- labxx
    }else{
        # pick three groups
        stopifnot(middleGroup %in% c(-1,0,1))
        
        lab <- rep(-1, length(x))
        labxx <- rep(0, length(xx))
        if(middleGroup == 0){
            labxx[fit3$dens[,1] > 99*fit3$dens[,3] &
                      fit3$dens[,1] > fit3$dens[,2]] <- -1
            labxx[fit3$dens[,3] > 99*fit3$dens[,1] &
                      fit3$dens[,3] > fit3$dens[,2]] <- 1
        }
        if(middleGroup == 1){
            labxx[fit3$dens[,1] > 99*fit3$dens[,3] &
                      fit3$dens[,1] > 99*fit3$dens[,2]] <- -1
            labxx[fit3$dens[,3] > 99*fit3$dens[,1] |
                      fit3$dens[,2] > 99*fit3$dens[,1]] <- 1
        }
        if(middleGroup == -1){
            labxx[fit3$dens[,1] > 99*fit3$dens[,3] |
                      fit3$dens[,2] > 99*fit3$dens[,3]] <- -1
            labxx[fit3$dens[,3] > 99*fit3$dens[,1] &
                      fit3$dens[,3] > 99*fit3$dens[,2]] <- 1
        }
        lab[x >= cut] <- labxx
    }
    
    return(list(label = lab, fit1 = fit1, fit2 = fit2, fit3 = fit3))
}

