library(ssadAssociation)

# ---------------

# I could simulate each spp separately, that would effectively make for more
# time steps for every iteration of the loop


nsite <- 20
nspp <- 50
allx <- replicate(nspp, numeric(nsite), simplify = FALSE)
allt <- numeric(nspp)
la <- 1.2
mu <- 1.45
rhom <- rbeta(nspp, 0.7, 1) * 1.5
rhol <- 0.15

nstep <- 8000

endtime <- 5000
xtot <- numeric(nstep)


# loop over steps
for(i in 1:nstep) {
    for(s in 1:nspp) {
        # only compute if this sp is not yet at target time
        if(allt[s] < endtime) {
            # print(s)
            # get a single sp
            x <- allx[[s]]

            # vector of rates:
            # birth = 1:nsite
            # death = nsite + (1:nsite)
            # imm locs = 2 * nsite + (1:nsite)
            # imm meta = 3 * nsite + 1
            r <- c(la * x, mu * x, rhol * x, rhom[s])

            # sample event
            e <- sample(1:length(r), 1, prob = r)

            # update based on event
            if(e <= nsite) { # birth
                iborn <- e
                x[iborn] <- x[iborn] + 1
            } else if(e <= 2 * nsite) { # death
                idead <- e - nsite
                x[idead] <- x[idead] - 1
            } else { # local imm or meta imm (effect is the same...for now)
                ireceive <- sample(nsite, 1)
                x[ireceive] <- x[ireceive] + 1
            }

            # log new pop state
            allx[[s]] <- x

            # update total time for this pop
            allt[s] <- allt[s] + rexp(1, sum(r))
        }
    }

    # track total individuals
    xtot[i] <- sum(unlist(allx))

    # stop if every spp reached target time
    if(all(allt >= endtime)) {
        print("all spp done")
        break
    }
}

if(any(allt < endtime)) {
    print("not all spp done")
    plot(allt - endtime)
}


plot(xtot, type = 'l')

# make site (row) by spp (col) matrix
sbys <- matrix(unlist(allx), nrow = nsite)
sbys


nbFit(sbys)

plusMinus(sbys)





# -----

# not sure what this is:

# I could simulate each spp sepparately, that would effectively make for more
# time steps for every iteration of the loop

metaN <- 10000
nsite <- 20
x <- numeric(nsite)
la <- 1.25
mu <- 1.45
rhom <- 1.25
rhol <- 0.2


# could make parameters like this:
#
# choose `mu` and `la + rhol` and `rhom` such that there's a stable pop after
# reasonable time and resulting in reasonable pop sizes
#
# assume a fixed ratio between `rhol` and `rhom` (relating to how proximate
# sites are and how distant metacomm is)
#
# fixed ratio of `rhol` and `rhom` will give value for `rhol` and `la`
#
# metacommunity imm should vary with spp based on metacomm abundance; under
# this scheme we should introduce `rhom0` which is per capita rate, such that
# `rhom0 * abund` = desired `rhom`


nstep <- 50000

latot <- mutot <- rholtot <- xtot <- numeric(nstep)


# loop over steps
for(i in 1:nstep) {
    # start with one species

    latot[i] <- sum(la * x)
    mutot[i] <- sum(mu * x)
    rholtot[i] <- sum(rhol * x)
    # vector of rates
    # birth = 1:nsite
    # death = nsite + (1:nsite)
    # imm locs = 2 * nsite + (1:nsite)
    # imm meta = 3 * nsite + 1
    r <- c(la * x, mu * x, rhol * x, rhom)

    # sample event
    e <- sample(1:length(r), 1, prob = r)

    # update based on event
    if(e <= nsite) { # birth
        iborn <- e
        x[iborn] <- x[iborn] + 1
    } else if(e <= 2 * nsite) { # death
        idead <- e - nsite
        x[idead] <- x[idead] - 1
    } else { # local imm or meta imm (effect is the same...for now)
        ireceive <- sample(nsite, 1)
        x[ireceive] <- x[ireceive] + 1
    }

    xtot[i] <- sum(x)
    # then combine with lapply

    # ultimately convert to C anyway
}

x
matrix(c(r[-length(r)], rep(r[length(r)], nsite)), ncol = 4)

plot(((rhom + rholtot + latot) / mutot)[1:10000], type = 'l', log = 'y')
abline(h = 1, col = 'red')

plot(xtot, type = 'l')




