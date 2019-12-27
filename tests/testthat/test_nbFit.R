context('nbFit works')

k <- 0.1
mu <- 5
set.seed(1)
r <- rnbinom(1000, k, mu = mu)
x <- cbind(r, r)
xfit <- nbFit(x)

test_that('nbFit returns correct dimension', {
    expect_equal(dim(xfit), c(2, 6))
})

test_that('nbFit returns correct parameter estimates', {
    expect_true(abs(xfit[1, 2] - k) < 0.05 * k)
    expect_true(abs(xfit[1, 3] - mu) < 0.05 * mu)
})
