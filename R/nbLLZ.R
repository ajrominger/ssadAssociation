#' @title Z-squared test statistic of the negative binomial distribution
#'
#' @description Exact test of the negative binomial
#'
#' @details Compares the observed log likelihood of the nagative binomial distribution
#' given the data to the sampling distribution of log likelihoods given the negative
#' binomial distribution is true.  This test statistic is Chi-square distributed with
#' d.f. = 1, meaning value larger than \code{qchisq(0.95, 1)} rejects the "null"
#' hypothesis that the data come from a negative binomial.
#'
#' @param x a numeric vector
#' @param size the size parameter
#' @param mu the mean parameter
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

nbLLZ <- function(x, size, mu) {
    likObs <- sum(dnbinom(x, size = size, mu = mu, log = TRUE))
    n <- length(x)

    # hypothetical distribution of probabilities
    p0 <- .p0(size, mu)

    # hypothetical mean and var
    m <- sum(p0 * exp(p0)) * n
    v <- sum((m/n - p0)^2 * exp(p0)) * n

    # z^2-value
    z <- ((likObs - m) / sqrt(v))^2

    return(as.numeric(z))
}


# @export

dlogLik <- function(x, size, mu, n, mod = c('nbinom', 'pois')) {
    mod <- match.arg(mod, c('nbinom', 'pois'))

    # hypothetical distribution of probabilities
    p0 <- .p0(size, mu, mod = mod)

    # hypothetical mean and var
    m <- sum(p0 * exp(p0)) * n
    v <- sum((m/n - p0)^2 * exp(p0)) * n

    return(dnorm(x, m, sqrt(v)))
}

# helper function to compute distribution of potential probabilities
.p0 <- function(size, mu, mod = 'nbinom') {
    n0 <- 0:10^5

    if(mod == 'nbinom') {
        p0 <- dnbinom(n0, size = size, mu = mu, log = TRUE)
    } else {
        p0 <- dpois(n0, lambda = mu, log = TRUE)
    }

    p0 <- p0[is.finite(p0)]

    if(exp(p0[length(p0)]) > .Machine$double.eps^0.75) {
        n0add <- (n0[length(p0)] + 1):10^6
        p0add <- dnbinom(n0add, size = size, mu = mu, log = TRUE)

        if(mod == 'nbinom') {
            p0add <- dnbinom(n0add, size = size, mu = mu, log = TRUE)
        } else {
            p0add <- dpois(n0add, lambda = mu, log = TRUE)
        }

        p0add <- p0add[is.finite(p0add)]

        p0 <- c(p0, p0add)
    }

    return(p0)
}
