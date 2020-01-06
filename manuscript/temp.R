N <- 4
nsite <- 30

klo <- replicate(1000,
                 max(schoener(
                     cbind(rnbinom(nsite, 0.1, N / nsite),
                           rnbinom(nsite, 0.1, mu = 1 / nsite)))))
klo <- replicate(1000,
                 max(schoener(
                     cbind(rnbinom(nsite, 100, N / nsite),
                           rnbinom(nsite, 100, mu = 1 / nsite)))))


plot(density(khi, na.rm = TRUE))
lines(density(klo, na.rm = TRUE), col = 'red')



image(kde2d(simNBPerm$nocc, simNBPerm$null_nocc), col = viridis(40))
points(simNBPerm[!duplicated(simNBPerm[, c('nocc', 'null_nocc')]), c('nocc', 'null_nocc')])

freq2d <- function(x, y, flim = NULL, trans = 'linear', ...) {
    f <- table(paste(x, y, sep = ','))

    xy <- lapply(strsplit(names(f), ','), as.numeric)
    xy <- do.call(rbind, xy)

    f <- as.numeric(f)

    # if(is.null(flim)) flim <- range(f)
    # plot(xy, pch = 21, bg = quantCol(f, viridis(50), trans = trans, xlim = flim), ...)

    o <- as.data.frame(cbind(xy, f))
    names(o) <- c('x', 'y', 'f')

    return(o)
}

freq2dPlot <- function(x, flim = NULL, trans = 'linear', ...) {
    if(is.null(flim)) flim <- range(x$f)
    plot(x$x, x$y, pch = 21,
         bg = quantCol(x$f, viridis(50), trans = trans, xlim = flim),
         ...)
}

foo1 <- freq2d(simNBPerm$nocc, simNBPerm$null_nocc, trans = 'log',
               xlim = range(rbind(simNBPerm, simPoPerm)[, c('nocc', 'null_nocc')]))

foo2 <- freq2d(simPoPerm$nocc, simPoPerm$null_nocc, trans = 'log',
               xlim = range(rbind(simNBPerm, simPoPerm)[, c('nocc', 'null_nocc')]))


freq2dPlot(foo1, flim = range(foo1$f, foo2$f), trans = 'log',
           xlim = c(1, nsite), ylim = c(1, nsite))
abline(0, 1, col = 'white', lwd = 4)
abline(0, 1, col = 'red', lwd = 1.5)

freq2dPlot(foo2, flim = range(foo1$f, foo2$f), trans = 'log')
abline(0, 1, col = 'white', lwd = 4)
abline(0, 1, col = 'red', lwd = 1.5)


lazyLoad('manuscript/RarePlusComMinus_cache/latex/ssad_perm_comp_053b819eb94a9bf6bde69c676716f3e4')

head(simNBPerm)
plot(simNBPerm$abund, simNBPerm$nocc, ylim = c(1, 30))
plot(simNBPerm$null_abund, simNBPerm$null_nocc, ylim = c(1, 30))

plot(simPoPerm$abund, simPoPerm$nocc, ylim = c(1, 30))
plot(simPoPerm$null_abund, simPoPerm$null_nocc, ylim = c(1, 30))


plot(simPoPerm$dAIC, simPoPerm$null_dAIC)
plot(simNBPerm$dAIC, simNBPerm$null_dAIC)


