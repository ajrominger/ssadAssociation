#' @title Fit the negative binomial
#'
#' @description Fits the negative binomial SSAD to a site x species matrix
#'
#' @details If likelihood optimization fails \code{NA} will be returned.
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
        fit <- try(MASS::fitdistr(x[, j], 'negative binomial'))

        if(inherits(fit, 'try-error')) {
            f <- c(size = NA, mu = NA, logLik = NA, poisLogLik = NA)
        } else {
            f <- c(fit$estimate, logLik = fit$loglik)
            ##### !!!!! use `mean(thisx[, j)` instead of `fit$estimate[2]`
            ##### !!!!! here and in the .Rmd file
            f <- c(f, poisLogLik = sum(dpois(x[, j], mean(x[, j]), log = TRUE)))
        }

        return(f)
    })

    theseStats <- t(rbind(spID = 1:ncol(x), theseStats))
    return(theseStats)
}
