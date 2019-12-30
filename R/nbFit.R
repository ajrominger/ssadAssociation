#' @title Fit the negative binomial
#'
#' @description Fits the negative binomial SSAD to a site x species matrix
#'
#' @details If likelihood optimization fails for a given species, \code{NA} will be
#' returned for all row elements corresponding to that species.
#'
#' @param x a site by species matrix
#'
#' @return A \code{matrix} with columns giving the species IDs, the two parameters of the
#' negative binomial distribution, the log likelihood of the negative binomial
#' distribution, and the log likelihood of the Poisson distribution.  Each row is one
#' species and the species ID refers to the column index of that species in the site x
#' species matrix \code{x}.
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

nbFit <- function(x) {
    theseStats <- sapply(1:ncol(x), function(j) {
        fit <- fitdistrNB(x[, j])
        f <- c(fit$estimate, logLik = fit$loglik)
        f <- c(f, poisLogLik = sum(dpois(x[, j], mean(x[, j]), log = TRUE)))
        f <- c(f, dAIC = as.numeric(2 * (2 - f['logLik'] - (1 - f['poisLogLik']))),
               z = nbLLZ(x[, j], f[1], f[2]),
               abund = sum(x[, j]), nocc = sum(x[, j] > 0))

        return(f)
    })

    theseStats <- t(rbind(spID = 1:ncol(x), theseStats))
    return(theseStats)
}
