library(pika)
library(vegan)

S0 <- 50
J <- 500
Jm <- 1e+07
mm <- 0.00001
ml <- 0.1
mu <- 1.2
la <- 1

# adding another site without increasing birth or imm rates is like decreasing
# the number of local individuals per site

nsite <- 3

nstep <- 50000

metaSAD <- rtnegb(S0, 200000, 1)
sum(metaSAD) - Jm

# site by spp matrix of local comms
X <- matrix(0, nrow = nsite, ncol = S0)

# vectors to track total richness and total N for all sites and spp
# (note: used to eval equilib)
totN <- numeric(nstep)
totS <- numeric(nstep)


# Z is normalization of rates to make probabilities below
Z <- mu * J + # total death rate
    la * J + # total dirth rate
    (nsite - 1) * J * ml + # total imm rate between local comms
    Jm * mm # total imm rate from metacomm

(pdead <- mu * J / Z)

(pbirth <- la * J / Z)

(pimmloc <- (nsite - 1) * J * ml / Z)

(pimmmet <- Jm * mm / Z)


for(i in 1:nstep) {
    # calculate death and birth rates per spp per site
    rdead <- X * mu
    rborn <- X * la

    # calculate local imm rate per species per site
    # (not this is rate of that spp in that site producing an immigrant
    # that will go somewhere else)
    rimmloc <- X * ml

    # calculate meta imm rate per spp
    rimmmet <- metaSAD * mm

    # combine into one vector
    r <- c(rdead, rborn, rimmloc, rimmmet)

    # choose event
    # 1:(S0 * nsite) = death
    # (S0 * nsite) + 1:(S * nsite) = birth
    # (2 * S0 * nsite) + 1:(S0 * nsite) = local imm
    # (3 * S0 * nsite) + 1:S0 = meta imm
    e <- sample(length(r), 1, prob = r)

    # determine outcome
    if(e <= S0 * nsite) { # death
        ijdead <- e
        X[ijdead] <- X[ijdead] - 1
    } else if(e <= 2 * S0 * nsite) { # birth
        ijborn <- e - (S0 * nsite)
        X[ijborn] <- X[ijborn] + 1
    } else if(e <= 3 * S0 * nsite) { # local imm
        ijimm <- e - (2 * S0 * nsite) # vector index
        ijimm <- arrayInd(ijimm, .dim = dim(X)) # array indeces

        # what site is it
        ifrom <- ijimm[1, 1]

        # what sp is it
        jsp <- ijimm[1, 2]

        # which site receives it
        ito <- sample((1:nrow(X))[-ifrom], 1)

        X[ito, jsp] <- X[ito, jsp] + 1
    } else { # meta imm
        iimm <- e - (3 * S0 * nsite)

        # which site receives it
        loc <- sample(nsite, 1)

        X[loc, iimm] <- X[loc, iimm] + 1
    }

    # update total N and S vectors
    totS[i] <- sum(colSums(X) > 0)
    totN[i] <- sum(X)
}

par(mfcol = c(3, 1), mar = c(2, 3, 0.5, 1), mgp = c(1.75, 0.5, 0))
plot(totS[ceiling(seq(1, nstep, length.out = 100))], type = 'l',
     ylab = 'S')
plot(totN[ceiling(seq(1, nstep, length.out = 100))], type = 'l',
     ylab = 'N')

par(mar = rep(0.5, 4), xaxs = 'i', yaxs = 'i')
Y <- X
Y[Y == 0] <- NA
image(t(Y), col = viridis::viridis(20), axes = FALSE, frame.plot = TRUE,
      panel.first = rect(-10, -10, S0 + 1, nsite + 1, col = 'gray'))

i <- 4
fitdistrNB(X[, i])
sum(dpois(X[, i], mean(X[, i]), log = TRUE))

mean(totN[(nstep - 1000):nstep])
