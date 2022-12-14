## 17/1/03. Written by Jens Schumacher. Please acknowledge as appropriate.

Xtree <- function(h)
    ## evaluate species branch matrix (sensu Petchey&Gaston) from a dendrogram
    ## tested for results of hclust and agnes
    ## hclust - hierarchical clustering 
    ## agnes - agglomerative clustering
    
    ## used components:
    ## merge - history of cluster merging
    ## height - actual heights at merging
    ## order - permutation to achieve nice output (needed only for agnes)
{

    species.names <- h$labels
    
    
    H1 <- matrix(0, length(h$order), 2 * length(h$order) - 2)
    l <- vector("numeric", 2 * length(h$order) - 2)
    for(i in 1:(length(h$order) - 1)) {
                                        # evaluate branch lengths
                                        #
        if(h$merge[i, 1] < 0) {
            l[2 * i - 1] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 1], 2 * i - 1] <- 1
        }
        else {
            l[2 * i - 1] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 1]]]
            H1[, 2 * i - 1] <- H1[, 2 * h$merge[i, 1] - 1] + H1[
                                                                , 2 * h$merge[i, 1]]
        }
        if(h$merge[i, 2] < 0) {
            l[2 * i] <- h$height[order(h$height)[i]]
            H1[ - h$merge[i, 2], 2 * i] <- 1
        }
        else {
            l[2 * i] <- h$height[order(h$height)[i]] - h$height[order(h$height)[h$merge[i, 2]]]
            H1[, 2 * i] <- H1[, 2 * h$merge[i, 2] - 1] + H1[, 2 *
                                                            h$merge[i, 2]]
        }
    }
    dimnames(H1) <- list(species.names,NULL)  
    list(h2.prime=l, H1=H1)
    ## l contains the length of all the tiny branches
    ## H1: each row represents one species, each column represents one branch
    ##     1 indicates that a branch is part of the pathway from species to top of the dendrogram
    ##     0 otherwise
}

