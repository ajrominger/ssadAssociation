#' @title Maximum likelihood of negative binomial distribution
#'
#' @description Finds the maximum likelihood estimates of the dispersion (k) and mean (mu)
#' parameters of the negative binomial distribution.
#'
#' @details This function lets \code{mu = mean(x)} and uses \code{optimize} to find
#' \code{k} with log likelihood function defined simply as \code{sum(dnbinom(x, k, mu))}
#'
#' @param x an integer vector
#'
#' @return A \code{lsit} with components \code{estiamte} giving the MLE and \code{loglik}
#' giving the value of the log likelihood function at the maximum. If
#' \code{var(x) < mean(x)} then the MLE does not exist and the funciton returns
#' \code{size = Inf}.
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

fitdistrNB <- function(x) {
    x <- x[!is.na(x)]

    fit <- optimise(.nbOptim, c(.Machine$double.eps, 1000), maximum = TRUE,
                    mu = mean(x), dat = x)

    names(fit) <- c('estimate', 'loglik')

    if(round(fit$estimate) == 1000) fit$estimate <- Inf

    fit$estimate <- c(size = fit$estimate, mu = mean(x))

    return(fit)
}


.nbOptim <- function(k, mu, dat) {
    sum(dnbinom(dat, size = k, mu = mu, log = TRUE))
}
