#' @title Duplicates \code{plusMinus} but using the independent swap algorithm
#' from \{picante\}
#'
#' @details All details are the same as \code{plusMinus} but this function uses
#' the independent swap algorithm from \{picante\} instead of the fixed-fixed
#' null algorithm. This function is only really intended to be used with
#' \code{indSwapTest}.
#'
#' @param x site by species matrix (sites as rows, species as columns, abundances in cells)
#' @param alpha the significance level
#' @param B number of replicates for the null model permutations
#'
#' @return an \code{ncol(x)} by \code{ncol(x)} matrix of species-species similarities
#'
#' @author Andy Rominger <ajrominger@@gmail.com>
#'
#' @export

plusMinusIndSwap <- function(x, alpha = 0.05, B = 999) {
    # metacommunity "abundnace"
    metaX <- colSums(x)

    # observed Schoener distance
    sd <- schoener(x)
    sd <- sd[lower.tri(sd)]

    # null Schoener distances
    nullsd <- lapply(1:B, function(i) {
        m <- randomizeMatrix(x, null.model = 'independentswap')
        d <- schoener(m)

        return(d[lower.tri(d)])
    })

    # combine null with observed
    nullsd <- cbind(sd, do.call(cbind, nullsd))

    # probabilities that the observed species-species distances are >= or <= the null
    ppos <- rowMeans(nullsd >= nullsd[, 1])
    pneg <- rowMeans(nullsd <= nullsd[, 1])

    # significantly positive and negative edges
    epos <- .signi(ppos, alpha)
    eneg <- .signi(pneg, alpha)

    # number of significantly positive and negative edges
    npos <- sum(epos)
    nneg <- sum(eneg)

    # edge list
    elist <- cbind(t(combn(1:ncol(x), 2)), epos, eneg)

    # correlation between centrality and abundance, and means
    corMeanPos <- .abundCenCorMean(elist[elist[, 3] == 1, 1:2, drop = FALSE], metaX)
    corMeanNeg <- .abundCenCorMean(elist[elist[, 4] == 1, 1:2, drop = FALSE], metaX)

    # number of species in each network
    vpos <- length(unique(as.vector(elist[elist[, 3] == 1, 1:2])))
    vneg <- length(unique(as.vector(elist[elist[, 4] == 1, 1:2])))

    return(list(all = c(v = ncol(x)),
                pos = c(n = npos, corMeanPos),
                neg = c(n = nneg, corMeanNeg)))
}
