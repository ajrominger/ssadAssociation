#' @title Simulate positive and negative interactions from purely random data
#'
#' @description Simulate spatial replicates of abundance data and apply the same
#' analytical pipeline to those data that would be applied to real data
#'
#' @details Note: any value passed to \code{ssadType} other than \code{'nbinom'} results
#' in a Poisson SSAD (i.e., there are only two options, negative binomial specified by
#' \code{'nbinom'} or Poisson specified by anything else)
#'
#' @param nsite number of sites to simulate
#' @param nspp number of species to simulate
#' @param mcCores number of cores to use in \code{parallel::mclapply}
#' @param ssadType string specifying SSAD shape (e.g. \code{'nbinom'})
#' @param sadfun function to generate random SAD sample
#' @param ssadfun function to generate random SSAD sample
#' @param nsim number of simulations to run
#'
#' @return an \code{ncol(x)} by \code{ncol(x)} matrix of species-species similarities
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export
#' @rdname simPlusMinus

simPlusMinus <- function(nsite, nspp, mcCores, ssadType, nsim) {
    return(NULL)
}


#' @export
#' @rdname simPlusMinus

simpleSim <- function(nsite, nspp, mcCores, sadfun, ssadfun, nsim) {
    # make SAD
    ii <- rep(1:nsim, each = nspp)
    X <- sadfun(nspp * nsim)
    X <- split(X, ii)

    o <- parallel::mclapply(X, mc.cores = mcCores, FUN = function(abund) {
        # calculate known quantities from abund
        J <- sum(abund)
        mu <- abund / nsite

        # loop over abundances and generate ssad
        mat <- sapply(1:length(abund), function(i) {
            return(ssadfun(nsite, mu[i]))
        })

        # browser()
        return(.simCleanup(mat))
    })

    o <- as.data.frame(do.call(rbind, o))
    o[is.na(o$pos.rho.rho) | is.na(o$neg.rho.rho), ] <- NA
    o <- o[!is.na(o$pos.n), ]

    return(o)
}

.simCleanup <- function(mat) {
    defaultNames <- c('pos.n', 'pos.rho.rho', 'pos.p', 'pos.m', 'pos.wm',
                      'neg.n', 'neg.rho.rho', 'neg.p', 'neg.m', 'neg.wm')
    mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

    if(any(dim(mat) < 10)) {
        o <- rep(NA, length(defaultNames))
    } else {
        o <- unlist(plusMinus(mat))
    }

    names(o) <- defaultNames

    return(o)
}

