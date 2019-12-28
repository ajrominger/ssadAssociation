context('fitdistrNB works like MASS::fitdistr')

test_that('fitdistrNB returns same results as MASS::fitdistr', {
    set.seed(1)
    x <- rnbinom(10000, 1, mu = 4)

    yoo <- unlist(fitdistrNB(x))
    bla <- unlist(MASS::fitdistr(x, 'negative binomial')[c('estimate', 'loglik')])

    d <- yoo - bla
    expect_true(all(abs(d) / ifelse(d < 0, bla, yoo) < .Machine$double.eps^0.25))
})
