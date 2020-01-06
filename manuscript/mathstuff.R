# foo <- replicate(100000, {
#     y <- sample(c(0, 0, 1))
#     r <- rpois(length(x), sum(x) / length(x))
#     z <- y
#     z[y == 1] <- 2
#
#     all(z == r)
# })

x <- c(5, 5, 0)

foo <- replicate(100000, x == rpois(length(x), mean(x)))

mean(foo)
prod(dpois(x, sum(x) / length(x)))


mean(foo) / prod(dpois(x, sum(x) / length(x)))



x <- c(5, 5, 5)

foo <- replicate(100000, {
    r <- rpois(length(x), mean(x))
    all(r == x)
})

mean(foo)

prod(dpois(x, mean(x)))



