context('schoener works')

test_that('RarePlusComMinus::schoener reproduces spaa::niche.overlap(x, method = "schoener")', {
    set.seed(1)
    x <- matrix(sample(10, 5 * 8, replace = TRUE), nrow = 5)

    mine <- schoener(x)
    pub <- as.matrix(spaa::niche.overlap(x, method = 'schoener'))

    expect_true(all.equal(as.vector(mine), as.vector(pub)))
})
