context('simulating plusMinus works without error')

nsim <- 2

test_that('simPlusMinus returns the correct dimension without error', {
    x <- simPlusMinus(nsitenspp = data.frame(nsite = rep(20, 2), nspp = rep(50, 2)),
                      sadStats = data.frame(mod = rep('fish', 2), par1 = rep(0.01, 2),
                                            par2 = rep(NA, 2)),
                      mcCores = 1, ssadType = 'nbinom',
                      kfun = function(nspp, abund) {return(0.5)}, nsim = 2)

    expect_equal(nsim, nrow(x))
    expect_equal(10, ncol(x))
})


test_that('simpleSim returns the correct dimension without error', {
    x <- simpleSim(nsite = 20, nspp = 100, mcCores = 1,
                   sadfun = function(n) {rfish(n, 0.01)},
                   ssadfun = function(nsite, mu) {rnbinom(nsite, 0.5, mu = mu)},
                   nsim = nsim)

    expect_equal(nsim, nrow(x))
    expect_equal(10, ncol(x))
})
