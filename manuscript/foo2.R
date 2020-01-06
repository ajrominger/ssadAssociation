b <- density(log(sadStats$par1[sadStats$mod == 'fish'], 10))
musig <- with(sadStats[sadStats$mod == 'plnorm', ], kde2d(par1, par2))
muk <- with(sadStats[sadStats$mod == 'tnegb', ], kde2d(par2, par1))

image(muk, col = viridis(20))

plot(b, type = 'n', xaxt = 'n')
segments(x0 = b$x[-1], x1 = b$x[-length(b$x)], y0 = b$y[-1], y1 = b$y[-length(b$y)],
         col = quantCol(b$y[-1], viridis(20)), lwd = 2)
axis(1, at = -3:0, labels = elabs(-3:0))

specialPoly <- function(x, y, col) {
    polygon(x, y, col = colAlpha(col, 0.4), border = NA)
    polygon(x, y, border = col)
}


hypRAD <- function(m, p, S) {
    x <- sad(model = m, par = p[!is.na(p)])

    r <- sad2Rank(x, S)

    return(r / sum(r))
}

S <- 100

allRAD <- mclapply(1:nrow(sadStats), mc.cores = nthrd, FUN = function(i) {
    hypRAD(sadStats$mod[i], as.numeric(sadStats[i, c('par1', 'par2')]), S)
})
allRAD <- do.call(rbind, allRAD)

radEnv <- lapply(unique(sadStats$mod), function(m) {
    apply(allRAD[sadStats$mod == m, ], 2, quantile, probs = c(0.025, 0.975))
})

plot(1, xlim = c(1, S), ylim = range(unlist(radEnv)), log = 'y', yaxt = 'n', type = 'n',
     xlab = 'Species rank', ylab = 'Relative abundance')
logAxis(2, expLab = TRUE)

polygon(c(1:S, S:1), c(radEnv[[1]][1, ], rev(radEnv[[1]][2, ])),
        col = hsv(0.56, 1, 0.8, 0.4))
polygon(c(1:S, S:1), c(radEnv[[2]][1, ], rev(radEnv[[2]][2, ])),
        col = hsv(0.05, 1, 1, 0.4))
polygon(c(1:S, S:1), c(radEnv[[3]][1, ], rev(radEnv[[3]][2, ])),
        col = hsv(0.75, 1, 0.8, 0.4))
polygon(c(1:S, S:1), c(radEnv[[1]][1, ], rev(radEnv[[1]][2, ])),
        border = hsv(0.56, 1, 0.8))
polygon(c(1:S, S:1), c(radEnv[[2]][1, ], rev(radEnv[[2]][2, ])),
        border = hsv(0.05, 1, 0.8))
polygon(c(1:S, S:1), c(radEnv[[3]][1, ], rev(radEnv[[3]][2, ])),
        border = hsv(0.75, 1, 0.8))

legend('topright',
       legend = c('Log-series', 'Poisson log-norm', 'zero-trunc. nbinom'),
       pch = 22, pt.cex = 2, pt.lwd = 1.5,
       col = hsv(c(0.56, 0.05, 0.75), 1, 0.8),
       pt.bg = hsv(c(0.56, 0.05, 0.75), 1, c(0.8, 1, 0.8), 0.4),
       bty = 'n')
