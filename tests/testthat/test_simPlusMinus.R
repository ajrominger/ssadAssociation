context('simPlusMinus works without error')

test_that('simPlusMinus returns the correct dimension without error', {
    nsim <- 3
    x <- simpleSim(nsite = 20, nspp = 100, mcCores = 1,
              sadfun = function(n) {rfish(n, 0.01)},
              ssadfun = function(nsite, mu) {rnbinom(nsite, 0.5, mu = mu)},
              nsim = nsim)

    expect_equal(nsim, nrow(x))
})
