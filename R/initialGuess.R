
# initial guess labels based on mixture model (truncated normal + normals)


.mixFit1 <- function(x){
    s0 <- sqrt(sum(x^2) / (length(x)-1))
    psingle <- dnorm(x, sd = s0, log = TRUE) + log(2)
    ll <- sum(psingle)
    return(list(pars = c(s0=s0), ll = ll, bic = -2*ll+1*log(length(x)),
                dens = matrix(psingle)))
}

.mixFit2 <- function(x, thresh = .01){
    # hidden factor, Z in [0,1], is prob observation X originated from left distn.
    # p1 := mean(Z)
    Z <- as.numeric(x <= max(x)*2/3)
    
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
        m2 <- weighted.mean(x, w = 1-Z) 
        s1 <- sqrt(sum(Z*(x^2)) / (sum(Z)-1))
        s2 <- weightedSd(x, 1-Z)
    }
    pleft <- dnorm(x, mean = 0, sd = s1) * p1 * 2
    pright <- dnorm(x, mean = m2, sd = s2) * (1-p1)
    S <- pleft + pright
    S[S==0] <- .Machine$double.eps
    ll <- sum(log(S))
    
    return(list(pars = c(p1=p1,s1=s1,m2=m2,s2=s2), ll = ll, bic = -2*ll+4*log(length(x)),
                dens = cbind(pleft, pright)))
}

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
        
        m2 <- weighted.mean(x, w = Z[,2])
        m3 <- weighted.mean(x, w = Z[,3])
        
        s1 <- sqrt(sum(Z[,1]*(x^2)) / (sum(Z[,1])-1))
        s2 <- weightedSd(x, Z[,2])
        s3 <- weightedSd(x, Z[,3])
    }
    pleft <- dnorm(x, mean = 0, sd = s1) * p[1] * 2
    pmid <- dnorm(x, mean = m2, sd = s2) * p[2]
    pright <- dnorm(x, mean = m3, sd = s3) * p[3]
    S <- pleft + pmid + pright
    S[S==0] <- .Machine$double.eps
    ll <- sum(log(S))
    
    return(list(pars = c(p1=p[1],p2=p[2],p3=p[3], m2=m2,m3=m3, s1=s1,s2=s2,s3=s3), 
                ll = ll, bic = -2*ll+7*log(length(x)),
                dens = cbind(pleft,pmid,pright)))
}



#' @param x The score (ie. debris score, doublet score, etc.) to be used for predicting each cell's label (eg. "doublet" vs. "cell"). 
#' @return A list with the following elements: \itemize{
#' \item{\code{label}} {A vector of the same length as \code{x} providing the labels (\code{-1} for cells, \code{1} for non-cells, \code{0} for uncertain).}
#' \item{\code{fit1}} {Summary of the 1-component (half Normal) model fit.}
#' \item{\code{fit2}} {Summary of the 2-component (half Normal + Normal) model fit.}
#' \item{\code{fit1}} {Summary of the 3-component (half Normal + 2 Normals) model fit.}}
initialGuess <- function(x){
    d <- density(x[which(x > min(x))]) 
    cut <- d$x[which.max(d$y)]
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
        labxx[fit2$dens[,1] > log(99)+fit2$dens[,2]] <- -1
        labxx[fit2$dens[,2] > log(99)+fit2$dens[,1]] <- 1
        lab[x >= cut] <- labxx
    }else{
        # pick three groups
        lab <- rep(-1, length(x))
        labxx <- rep(0, length(xx))
        labxx[fit3$dens[,1] > log(99)+fit3$dens[,2]] <- -1
        labxx[fit3$dens[,3] > fit3$dens[,1] &
                  fit3$dens[,3] > fit3$dens[,2]] <- 1
        lab[x >= cut] <- labxx
    }
    
    return(list(label = lab, fit1 = fit1, fit2 = fit2, fit3 = fit3))
}

plotInitialGuess <- function(x, IG = NULL, fit = NULL){
    
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
    
    layout(matrix(1:2, nrow=1))
    hist(x, breaks=100, probability = TRUE, main = 'Full Score')
    abline(v = cut)
    if(fit == 1){
        curve(p0*2*dnorm(x-cut,sd=IG$fit1$pars[1]), from = cut, to = max(x)+1, col = 4, lwd=2, add = TRUE)
    }
    if(fit == 2){
        p1 <- IG$fit2$pars[1]; s1 <- IG$fit2$pars[2]; m2 <- IG$fit2$pars[3]; s2 <- IG$fit2$pars[4]
        curve(p0*2*p1*dnorm(x-cut,sd=s1), from = cut, to = max(x)+1, col = 4, lwd=2, add = TRUE)
        curve(p0*(1-p1)*dnorm(x-cut,mean=m2,sd=s2), from = cut, to = max(x)+1, col = 2, lwd=2, add = TRUE)
    }
    if(fit == 3){
        p1 <- IG$fit3$pars[1]; p2 <- IG$fit3$pars[2]; p3 <- IG$fit3$pars[3]
        m2 <- IG$fit3$pars[4]; m3 <- IG$fit3$pars[5]
        s1 <- IG$fit3$pars[6]; s2 <- IG$fit3$pars[7]; s3 <- IG$fit3$pars[8]
        curve(p0*2*p1*dnorm(x-cut,sd=s1), from = cut, to = max(x)+1, col = 4, lwd=2, add = TRUE)
        curve(p0*p2*dnorm(x-cut,mean=m2,sd=s2), from = cut, to = max(x)+1, col = 1, lwd=2, add = TRUE)
        curve(p0*p3*dnorm(x-cut,mean=m3,sd=s3), from = cut, to = max(x)+1, col = 2, lwd=2, add = TRUE)
    }
    
    hist(xx, breaks=100, probability = TRUE, main = 'Truncated Score')
    if(fit == 1){
        curve(2*dnorm(x,sd=IG$fit1$pars[1]), from = 0, to = max(x)+1, col = 4, lwd=2, add = TRUE)
    }
    if(fit == 2){
        p1 <- IG$fit2$pars[1]; s1 <- IG$fit2$pars[2]; m2 <- IG$fit2$pars[3]; s2 <- IG$fit2$pars[4]
        curve(2*p1*dnorm(x,sd=s1), from = 0, to = max(xx)+1, col = 4, lwd=2, add = TRUE)
        curve((1-p1)*dnorm(x,mean=m2,sd=s2), from = 0, to = max(xx)+1, col = 2, lwd=2, add = TRUE)
    }
    if(fit == 3){
        p1 <- IG$fit3$pars[1]; p2 <- IG$fit3$pars[2]; p3 <- IG$fit3$pars[3]
        m2 <- IG$fit3$pars[4]; m3 <- IG$fit3$pars[5]
        s1 <- IG$fit3$pars[6]; s2 <- IG$fit3$pars[7]; s3 <- IG$fit3$pars[8]
        curve(2*p1*dnorm(x,sd=s1), from = 0, to = max(xx)+1, col = 4, lwd=2, add = TRUE)
        curve(p2*dnorm(x,mean=m2,sd=s2), from = 0, to = max(xx)+1, col = 1, lwd=2, add = TRUE)
        curve(p3*dnorm(x,mean=m3,sd=s3), from = 0, to = max(xx)+1, col = 2, lwd=2, add = TRUE)
    }
    layout(1)
}

