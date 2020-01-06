lazyLoad('manuscript/RarePlusComMinus_cache/latex/obs_plusMinus_comp_347ada9e8b34620896b2e0ae7240c95a')
lazyLoad('manuscript/RarePlusComMinus_cache/latex/sim_plusMinus_comp_b06e34a53f96a43ccd3e563d692743fc')

hist(commStats$pos.m - commStats$neg.m, breaks = seq(-0.3, 0.3, by = 0.1/3))
hist(commStats$pos.wm - commStats$neg.wm, breaks = seq(-0.3, 0.3, by = 0.1/3))

specialHist(commStats$pos.m, breaks = seq(0, 0.3, by = 0.1/6), col = hsv(0.6, alpha = 0.5))
specialHist(commStats$neg.m, breaks = seq(0, 0.3, by = 0.1/6), col = hsv(0, alpha = 0.5), add = TRUE)

specialHist(simPMData$pos.m, breaks = seq(0, 0.3, by = 0.1/6), col = hsv(0.6, alpha = 0.5))
specialHist(simPMData$neg.m, breaks = seq(0, 0.3, by = 0.1/6), col = hsv(0, alpha = 0.5), add = TRUE)

lazyLoad('manuscript/RarePlusComMinus_cache/latex/sad_comp_e72ffcd46a05da0d56045916fe8fecde')

with(sadStats[sadStats$mod == 'fish', ],
     plot(nsite, par1))
boxplot(nspp ~ mod, data = sadStats)


lazyLoad('manuscript/RarePlusComMinus_cache/latex/sim_simple_pois_plusMinus_comp_9aec93b034a294a1d0656ce5cb025984')
lazyLoad('manuscript/RarePlusComMinus_cache/latex/sim_simple_nb_plusMinus_comp_3d6d4bd446c13a2e831449802d26b715')
lazyLoad('manuscript/RarePlusComMinus_cache/latex/sim_pois_plusMinus_comp_9a96fbe13b714de0ae17332fbc684f1d')

plot(sort(simPMSimpNB$pos.n), ylim = c(0, 250))
points(sort(simPMSimpPo$pos.n), col = 'red')

plot(sort(simPMSimpNB$neg.n), ylim = c(0, 2500))
points(sort(simPMSimpPo$neg.n), col = 'red')




schoener(cbind(c(1, 1, 0, 0), c(1, 1, 0, 0)))


schoenerSim <- function(nsite, nspp, mcCores, sadfun, ssadfun, nsim) {
    # make SAD
    ii <- rep(1:nsim, each = nspp)
    X <- sadfun(nspp * nsim)
    X <- split(X, ii)

    o <- parallel::mclapply(X, mc.cores = mcCores, FUN = function(abund) {
    # o <- lapply(X, function(abund) {
        # calculate known quantities from abund
        # J <- sum(abund)
        mu <- abund / nsite

        # loop over abundances and generate ssad
        mat <- .makeMat(nsite, nspp, ssadfun, mu = mu)
        mat <- mat[rowSums(mat) > 0, colSums(mat) > 0]
        abund <- colSums(mat) / sum(mat)

        # run null model
        nulls <- simulate(nullmodel(mat, 'r2dtable'), nsim = 1)

        # calculate schoener
        schoObs <- schoener(mat)
        schoPerm  <- schoener(nulls[, , 1])
        ij <- t(combn(1:ncol(mat), 2))

        o <- cbind(abund[ij[1, ]], abund[ij[2, ]],
                   schoObs[lower.tri(schoObs)], schoPerm[lower.tri(schoPerm)])
        colnames(o) <- c('abund1', 'abund2', 'schoObs', 'schoPerm')

        return(o)
    })

    return(as.data.frame(do.call(rbind, o)))
}

schoSim <- schoenerSim(30, 50, nthrd, function(n) {rfish(n, 0.01)},
                       function(n, mu) {rnbinom(n, 0.1, mu = mu)}, nsim = 50)

plot((schoSim$abund1 + schoSim$abund2) / 2, schoSim$schoObs - schoSim$schoPerm)



