library(pika)

S <- 150
b <- 0.001
nsite <- 25

x <- rfish(S, b)
y <- rnbinom(nsite * length(x), size = 0.1, mu = rep(x, nsite) / nsite)
y <- t(matrix(y, ncol = nsite))
y <- y[rowSums(y) > 0, ]

pos.neg.abun.r2d(y, 'foo', 0.05, plots = FALSE)

z <- y / rowSums(y)
as.matrix(dist(z[1, ]))

dim(outer(z, z, '-'))

plot(sort(x), sort(z), xlim = range(x, z), ylim = range(x, z))
abline(0, 1, col = 'red')


plot(dnbinom(0:100, size = 0.8, mu = 10), log = 'y')
