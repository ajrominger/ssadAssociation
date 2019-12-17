# library(pika)
#
# S <- 4
# b <- 0.001
# nsite <- 3
#
# x <- rfish(S, b)
# y <- rnbinom(nsite * length(x), size = 100, mu = rep(x, nsite) / nsite)
# y <- t(matrix(y, ncol = nsite))
#
#
# z <- y / rowSums(y)
# as.matrix(dist(z[1, ]))
#
# dim(outer(z, z, '-'))
#
# plot(sort(x), sort(z), xlim = range(x, z), ylim = range(x, z))
# abline(0, 1, col = 'red')
#
#
# plot(dnbinom(0:100, size = 0.8, mu = 10), log = 'y')
