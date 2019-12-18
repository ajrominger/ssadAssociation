context('plusMinus reproduces results of published function')

library(spaa)
library(igraph)
library(vegan)

# function from XYZ, except that:
# 1) packages are loaded outside the function, not inside
# 2) results are not written out (the line `save(resFin,file...` is commented out)

pos.neg.abun.r2d <- function(matb,name,alpha,plots){

    #INPUT:
    # matb = abundance matrix with species in colunms
    # name = name to save results
    # alpha = nominal error to detected significant associations
    # plots = T or F for ploting networks with community information
    #
    #STEPS:
    #a)Network generations
    #	1:	Using the Schoener's index to measure similarity in the distribution of abundances across sample units (quadrats, plots, pitfalls etc. ) for each pair of species.
    #	2:	Assessing significance by a fixed-fixed null model (?r2dtable? in Vegan package)
    #	3:	Selecting significantly aggregated species pairs to depict a link in the positive network.
    #	4:	Selecting significantly segregated species pairs to depict a link in the negative network
    #	5:	Estimate network properties

    #b) calculations:
    #	1: "num" = number of links
    #	2: "num_nodos" = number of nodes
    #	3: "con" = connectivity
    #	4: "mod_val" = modularity (best modular partition over 100 runs of Louvain)
    #	5: "abu_rho" = Spearman's rho coefficient between species abundances and degree
    #	6: "abu_p" = p-value of Spearman's rho coefficient between species abundances and degree
    #	7: "w.mean.abu" = mean species abundances weighted by species degree



    #supplementary function to asses significance
    signi<-function(x,alpha){y <- x; y[x>alpha]<-0; y[x<=alpha] <- 1; return(y)}


    matb<-matb[rowSums(matb)>0,colSums(matb)>0]

    abu <- colSums(matb)
    #species abundance

    Sd<-niche.overlap(matb,method="schoener")
    # Schoener overlap

    nul_mats <- list()
    for(i in 1:999){
        m <- nullmodel(matb, "r2dtable")
        nul_mats[[i]]<-niche.overlap(simulate(m, nsim=1)[,,1],method="schoener")
    }
    nul_res <- do.call(rbind,nul_mats)
    # null Schoener overlap
    browser()

    nul_res<-rbind(Sd,nul_res)
    pos_p <- apply(nul_res,2,function(x){sum(x>=x[1])/length(x)})
    neg_p <- apply(nul_res,2,function(x){sum(x<=x[1])/length(x)})


    pos_casos <- signi(pos_p,alpha)
    neg_casos <- signi(neg_p,alpha)
    # significant pairs

    num_pos <- sum(pos_casos)
    num_neg <- sum(neg_casos)

    verlist <- data.frame(t(combn(colnames(matb),2)),pos_casos,neg_casos)
    pos <- verlist[verlist[,3]==1,1:2]
    neg <- verlist[verlist[,4]==1,1:2]
    red_pos <- graph_from_data_frame(pos, directed=F)
    red_neg <- graph_from_data_frame(neg, directed=F)
    # positives and negatives networks

    if(num_pos>0){
        lmod<-list()
        for(i in 1:100){
            lmod[[i]] <- cluster_louvain(red_pos, weights = NULL)
        }
        mod_pos <- unlist(lapply(lmod,function(x)max(x$modularity)))
        lmod_pos<-lmod[[sample(which(mod_pos==max(mod_pos)),1)]]
        mod_val_pos<-max(lmod_pos$modularity)
    } else {mod_val_pos <- "NA"}
    # modularity positive

    if(num_neg>0){
        lmod<-list()
        for(i in 1:100){
            lmod[[i]] <- cluster_louvain(red_neg, weights = NULL)
        }
        mod_neg <- unlist(lapply(lmod,function(x)max(x$modularity)))
        lmod_neg<-lmod[[sample(which(mod_neg==max(mod_neg)),1)]]
        mod_val_neg<-max(lmod_neg$modularity)
    } else {mod_val_neg <- "NA"}
    # modularity neg


    num_nodos_pos <- length(unique(c(as.character(pos[,1]),as.character(pos[,2]))))
    num_nodos_neg <- length(unique(c(as.character(neg[,1]),as.character(neg[,2]))))
    #number of nodes

    con_pos <- num_pos/choose(num_nodos_pos,2)
    con_neg <- num_neg/choose(num_nodos_neg,2)
    # network connectivity

    cen_pos <-table(c(as.character(pos[,1]),as.character(pos[,2])))
    cen_neg <-table(c(as.character(neg[,1]),as.character(neg[,2])))
    abu_cor_pos <- abu[names(abu)%in%names(cen_pos)]
    abu_cor_pos <- abu_cor_pos[match(names(cen_pos),names(abu_cor_pos))]
    abu_cor_neg <- abu[names(abu)%in%names(cen_neg)]
    abu_cor_neg <- abu_cor_neg[match(names(cen_neg),names(abu_cor_neg))]
    cor_abu_cen_pos <- cor.test(abu_cor_pos,cen_pos,method="spearman")
    cor_abu_cen_neg <- cor.test(abu_cor_neg,cen_neg,method="spearman")
    abu_rho_pos <- cor_abu_cen_pos$estimate
    abu_rho_neg <- cor_abu_cen_neg$estimate
    abu_p_rho_pos <- cor_abu_cen_pos$p.value
    abu_p_rho_neg <- cor_abu_cen_neg$p.value
    #abudance centrality correlations

    pos.mean.abu<-mean (abu_cor_pos/sum(matb))
    neg.mean.abu<-mean (abu_cor_neg/sum(matb))
    w.pos.mean.abu<-weighted.mean(abu_cor_pos/sum(matb),cen_pos)
    w.neg.mean.abu <-weighted.mean(abu_cor_neg/sum(matb),cen_neg)
    #weighted mean relative abundances


    resFin<-cbind(name,num_pos,num_neg,num_nodos_pos,num_nodos_neg,con_pos,con_neg,mod_val_pos,mod_val_neg,
                  abu_rho_pos,abu_rho_neg,abu_p_rho_pos ,abu_p_rho_neg,w.pos.mean.abu, w.neg.mean.abu)
    # final results

    # save(resFin,file= paste(name,"_alpha_",alpha,".Rdata",sep=""))
    # save results


    if(plots){
        if(num_pos>1 & num_neg>1){
            pdf(paste(name,"_alpha_",alpha,".pdf",sep=""),width=12,height=6)
            par(mfrow=c(1,2))
            V(red_pos)$label.cex = 0.7
            V(red_neg)$label.cex = 0.7
            l2 <- layout.fruchterman.reingold(red_pos)
            plot(red_pos,vertex.color=lmod_pos$membership,layout=l2)
            text(-1,1.3,paste(name,"   Positive","Q =",round(mod_val_pos,2)))
            l2 <- layout.fruchterman.reingold(red_neg)
            plot(red_neg,vertex.color=lmod_neg$membership,layout=l2)
            text(-1,1.3,paste("   Negative","Q =",round(mod_val_neg,2)))
            dev.off()
        }
        #save & plot networks + community structure
    }
    return(resFin)
}


# test_that('plusMinus returns same results as published version', {
#     pos.neg.abun.r2d(matb = abun.mat[[1]], name = 'foo', alpha = 0.05, plots = FALSE)
#     set.seed(1)
#     x <- matrix(sample(10, 5 * 8, replace = TRUE), nrow = 5)
#
#     mine <- schoener(x)
#     pub <- as.matrix(spaa::niche.overlap(x, method = 'schoener'))
#
#     expect_true(all.equal(as.vector(mine), as.vector(pub)))
# })
