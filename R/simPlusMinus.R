#' @title Simulate positive and negative interactions from purely random data
#'
#' @description Simulate spatial replicates of abundance data and apply the same
#' analytical pipeline to those data that would be applied to real data
#'
#' @details \code{simPlusMinus} draws random SAD and SSAD shapes from the raw data
#' and uses these to simulate more data and calculate network statistics on those
#' simulated data. \code{simpleSim} assumes one SAD and one SSAD and simulates data from
#' those, again calculated network statistics.
#'
#' Note: any value passed to \code{ssadType} other than \code{'nbinom'} results
#' in a Poisson SSAD (i.e., there are only two options, negative binomial specified by
#' \code{'nbinom'} or Poisson specified by anything else)
#'
#' @param sadStats a \code{data.frame} with columns \code{mod}, \code{par1}, \code{par2}
#' @param nsite number of sites to simulate
#' @param nspp number of species to simulate
#' @param mcCores number of cores to use in \code{parallel::mclapply}
#' @param ssadType string specifying SSAD shape (e.g. \code{'nbinom'})
#' @param sadfun function to generate random SAD sample
#' @param ssadfun function to generate random SSAD sample
#' @param kfun function to relate k parameter of the SSAD to abundance
#' @param nsim number of simulations to run
#'
#' @return a \code{data.frame} with \code{<= nsim} rows (some simulations may be
#' thrown out if they do not meet data filtering requirements), and columns corresponding
#' to summary statistics about the positive and negative network characteristics
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export
#' @rdname simPlusMinus

simPlusMinus <- function(sadStats, mcCores,
                         ssadType = 'nbinom', kfun, nsim) {
    # indeces for SAD data and nsite, nspp data
    iiSAD <- sample(nrow(sadStats), nsim, replace = TRUE)
    # jjNN <- sample(nrow(nsitenspp), nsim, replace = TRUE)

    # loop over replicates and make communities, then calculate networks
    o <- parallel::mclapply(1:nsim, mc.cores = mcCores, FUN = function(i) {
        j <- iiSAD[i]
        nsite <- sadStats$nsite[j]
        nspp <- sadStats$nspp[j]

        # make SAD
        # iSAD <- iiSAD[i]
        rfun <- get(paste0('r', sadStats$mod[j]))
        pars <- as.numeric(sadStats[j, 2:3])
        pars <- pars[!is.na(pars)]
        abund <- do.call(rfun, c(list(nspp), as.list(pars)))

        # calculate known quantities from abund
        J <- sum(abund)
        mu <- abund / nsite

        # loop over abundances and generate ssad
        if(ssadType == 'nbinom') {
            # calcualte k (size param)
            k <- kfun(nspp, abund)

            mat <- .makeMat(nsite, nspp, rnbinom, size = k, mu = mu)
        } else {
            mat <- .makeMat(nsite, nspp, rpois, lambda = mu)
        }

        return(.simCleanup(mat))
    })

    return(.outCleanup(o))
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
        mat <- .makeMat(nsite, nspp, ssadfun, mu = mu)

        return(.simCleanup(mat))
    })

    return(.outCleanup(o))
}


.makeMat <- function(nsite, nspp, rfun, ...) {
    matrix(rfun(nsite * nspp, ...), nrow = nsite, byrow = TRUE)
}


.simCleanup <- function(mat) {
    defaultNames <- c('all.v',
                      'pos.n', 'pos.v', 'pos.rho.rho', 'pos.p', 'pos.m', 'pos.wm',
                      'neg.n', 'neg.v', 'neg.rho.rho', 'neg.p', 'neg.m', 'neg.wm')
    mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

    if(any(dim(mat) < 10)) {
        o <- rep(NA, length(defaultNames))
    } else {
        o <- unlist(plusMinus(mat))
    }

    names(o) <- defaultNames

    return(o)
}

.outCleanup <- function(o) {
    o <- as.data.frame(do.call(rbind, o))
    o[is.na(o$pos.rho.rho) | is.na(o$neg.rho.rho), ] <- NA
    o <- o[!is.na(o$pos.n), ]

    return(o)
}
