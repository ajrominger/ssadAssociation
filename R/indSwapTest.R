#' @title Explore the effect of using the independent swap algorithm on inference
#' of connection between abundance and positive or negative association networks
#'
#' @description This is largely a copy of the function \code{simPlusMinus} but
#' using the independent swap algorithm from \{picante\}
#'
#'
#' @param sadStats a \code{data.frame} with columns \code{mod}, \code{par1}, \code{par2}
#' @param nsite number of sites to simulate
#' @param nspp number of species to simulate
#' @param mcCores number of cores to use in \code{parallel::mclapply}
#' @param ssadType string specifying SSAD shape (e.g. \code{'nbinom'})
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

indSwapTest <- function(sadStats, mcCores,
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

        return(.simCleanupIndSwap(mat))
    })

    return(.outCleanup(o))
}


.simCleanupIndSwap <- function(mat) {
    defaultNames <- c('all.v',
                      'pos.n', 'pos.v', 'pos.rho.rho', 'pos.p', 'pos.m', 'pos.wm',
                      'neg.n', 'neg.v', 'neg.rho.rho', 'neg.p', 'neg.m', 'neg.wm')
    mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

    if(any(dim(mat) < 10)) {
        o <- rep(NA, length(defaultNames))
    } else {
        o <- unlist(plusMinusIndSwap(mat))
    }

    names(o) <- defaultNames

    return(o)
}
