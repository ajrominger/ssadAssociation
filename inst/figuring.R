#' @title Positive and negative associations
#'
#' @description Calculate positive and negative associations based on Schoener similarity
#'
#' @details Faster implementation of the funciton \code{pos.neg.abun.r2d} from ((citation))
#' that also only returns information relevant to testing the association of abundance and
#' different types of connectivity (plus and minus).
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

plusMinus <- function(x, alpha = 0.05, B = 999) {
    x <- x[rowSums(x) > 0, colSums(x) > 0]

    # metacommunity "abundnace"
    metaX <- colSums(x)

    # observed Schoener distance
    sd <- as.vector(schoener(x))

    # null Schoener distances
    nulls <- simulate(nullmodel(x, 'r2dtable'), nsim = B)
    nullsd <- lapply(1:dim(nulls)[3], function(i) {
        return(as.vector(schoener(nulls[, , i])))
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
    elist <- data.frame(expand.grid(1:ncol(x), 1:ncol(x)), epos, eneg)
    # pos <- verlist[verlist[, 3] == 1, 1:2]
    # neg <- verlist[verlist[, 4] == 1, 1:2]
    # gpos <- graph_from_data_frame(verlist[verlist[, 3] == 1, 1:2],
    #                                  directed = FALSE)
    # gneg <- graph_from_data_frame(verlist[verlist[, 4] == 1, 1:2],
    #                                  directed = FALSE)

    # correlation between centrality and abundance, and means
    corMeanPos <- .abundCenCorMean(elist[elist[, 3] == 1, 1:2], metaX)
    corMeanNeg <- .abundCenCorMean(elist[elist[, 4] == 1, 1:2], metaX)

    return(list(pos = c(n = npos, corMeanPos),
                neg = c(n = nneg, corMeanPos)))
}


# function to extract probabilities <= alpha
.signi <- function(x, alpha) {
    return(as.integer(x <= alpha))
}


# function to calculate correlation between abundance and centrality
# and mean abundances, unweighted and weighted by centrality
.abundCenCorMean <- function(elist, abund) {
    cen <- table(c(elist[, 1], elist[, 2]))
    x <- abund[as.integer(names(cen))]

    # correlation output
    o <- cor.test(x, cen, method = 'spearman')

    # relative abundance
    relx <- x / sum(abund)

    # return: corr coeff, corr p-val, mean abund, mean abund weighted by centrality
    return(list(rho = o$estimate, p = o$p.value,
                m = mean(relx), wm = weighted.mean(relx, cen)))
}
