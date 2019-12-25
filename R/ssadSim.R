#' @title Simulate permuted SSADs
#'
#' @description Simulate many SSADs and permute them to understand the effect of
#' permutation on the shape of the SSAD
#'
#' @details Note, \code{sadfun} must be a function that takes a single argument \code{n}
#' (number of species), e.g., \code{function(n) pika::rfish(n, b = 0.1)}. Similarely,
#' \code{ssadfun} must take both \code{n} (the number of sites) and also \code{mu} the
#' mean abundance, e.g. \code{function(n) rnbinom(n, size = 0.1, mu = mu)}.
#'
#' @param nsite number of sites to simulate
#' @param nspp number of species to simulate
#' @param mcCores number of cores to use in \code{parallel::mclapply}
#' @param sadfun function to generate random SAD sample
#' @param ssadfun function to generate random SSAD sample
#' @param nsim number of simulations to run
#' @param B number of replicates for the null model permutations
#'
#' @return a \code{data.frame} with columns giving the z-score of the two negative
#' binomial parameters, and the difference in AIC values of the negative binomial versus
#' the Poisson
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

ssadSim <- function(nsite, nspp, mcCores, sadfun, ssadfun, nsim, B = 999) {
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
        mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]

        # run null model
        nulls <- simulate(nullmodel(mat, 'r2dtable'), nsim = B)

        # calculate ssad fit info
        nullfit <- lapply(1:dim(nulls)[3], function(i) {
            return(nbFit(nulls[, , i]))
        })

        nullfit <- do.call(rbind, nullfit)

        nullMSD <- lapply(1:ncol(mat), function(i) {
            nulldat <- nullfit[nullfit[, 1] == i, ]
            dAIC <- 2 * ((2 - nulldat[, 4]) - (1 - nulldat[, 5]))

            c(sizeM = mean(nulldat[, 2]), sizeSD = sd(nulldat[, 2]),
              muM = mean(nulldat[, 3]), muSD = sd(nulldat[, 3]),
              dAICM = mean(dAIC), dAICSD = sd(dAIC))
        })
        nullMSD <- do.call(rbind, nullMSD)

        obs <- nbFit(mat)
        z <- cbind(size = .z(obs[, 2], nullMSD[, 1], nullMSD[, 2]),
                   mu = .z(obs[, 3], nullMSD[, 3], nullMSD[, 4]),
                   dAIC = .z(2 * ((2 - obs[, 4]) - (1 - obs[, 5])),
                             nullMSD[, 5], nullMSD[, 6]))
    })

    return(as.data.frame(do.call(rbind, o)))
}

.z <- function(x, m, s) {
    s[s < .Machine$double.eps] <- 1

    (x - m) / s
}
