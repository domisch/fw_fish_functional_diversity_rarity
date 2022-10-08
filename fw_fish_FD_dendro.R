fw_fish_FD_dendro <- function (S, A, w = NA, Distance.method = "gower", ord = c("podani", "metric"), 
                          Cluster.method = c(ward = "ward", single = "single", 
                                             complete = "complete", UPGMA = "average", UPGMC = "centroid", 
                                             WPGMC = "median", WPGMA = "mcquitty"), stand.x = TRUE, stand.FD = FALSE, 
                          Weigthedby = c("abundance", "biomasCarabids", "biomasBees", 
                                         "biomassValue"), biomassValue = NA) 
{
  require(FD)
  require(cluster)
  require(vegan)
  Out <- data.frame(comm = rep(NA, nrow(A)), n_sp = rep(NA, 
                                                        nrow(A)), n_tr = rep(NA, nrow(A)), FDpg = rep(NA, nrow(A)), 
                    FDw = rep(NA, nrow(A)), FDwcomm = rep(NA, nrow(A)), 
                    qual.FD = rep(NA, nrow(A)))
  Out$comm <- as.numeric(rownames(A))
  Out$n_tr <- ncol(S)
  Arich <- as.matrix(A)
  Arich[which(Arich > 0)] <- 1
  Out$n_sp <- rowSums(Arich, na.rm = TRUE)
  if (is.na(w)[1]) {
    w <- rep(1, ncol(S))
  }
  if (Distance.method == "gower") {
    D <- gowdis(S, w = w, ord = ord)
  }  else {
    if (stand.x == TRUE) {
      S2 <- scale(S, center = TRUE, scale = TRUE)
      D <- dist(S2, method = Distance.method)
    }  else {
      D <- dist(S, method = Distance.method)
    }
  }
  tree <<- hclust(D, method = Cluster.method)
  # plot(tree)
  xtree <- Xtree(tree)
  c_distance <- cor(D, cophenetic(tree))
  Out[, 7] <- rep(c_distance, nrow(Out))
  AA <- A
  if (Weigthedby != "abundance") {
    if (Weigthedby == "biomasCarabids") {
      biomassValue2 <- Jelaska(biomassValue)
    }
    if (Weigthedby == "biomasBees") {
      biomassValue2 <- Cane(biomassValue)
    }  else {
      biomassValue2 <- biomassValue
    }
    if (is.vector(biomassValue2)) {
      for (j in 1:ncol(A)) AA[, j] <- A[, j] * biomassValue2[j]
    } else {
      AA <- AA * biomassValue2
    }
  }
  AFw <- AA
  AFw <- as.data.frame(t(apply(AA, 1, function(x) {x/max(x)} )))
  # for (k in 1:nrow(AA)) {
  #   AFw[k, ] <- AA[k, ]/max(AA[k, ])
  # }
  AFcomm <- AA
  tmp <- max(AA)
  AFcomm <- as.data.frame(t(apply(AA, 1, function(x) {x/tmp} )))
  # for (k in 1:nrow(AA)) {
  #   AFcomm[k, ] <- AA[k, ]/max(AA)
  # }
  
  ###---- Use foreach -----#
  ### Sends only one row 
  Out <-  foreach(i=1:nrow(A), fw_fish_row=mylist, .inorder = T, .combine = 'rbind', .export=c("xtree"), 
                  .verbose=F, .errorhandling=c('stop'),
                  .packages=c("fundiv", "FD", "cluster", "vegan")) %dopar% {
                    
                #   for (i in 1:nrow(A)) {                        # original
                #   species_in_C <- ifelse(A[i, ] > 0, 1, 0)      # original
                    species_in_C <- ifelse(fw_fish_row > 0, 1, 0)
                    select_xtree <- xtree$H1[which(species_in_C > 0), ]
                    if (is.array(select_xtree) == TRUE) {
                      i.primeC <- ifelse(colSums(select_xtree) > 0, 1,0)
                    }  else {
                      i.primeC <- select_xtree
                    }
                    Out[i, 4] <- sum(i.primeC * xtree$h2.prime)
                    xtree.weigths <- xtree$H1
                    for (k in 1:nrow(S)) {
                      xtree.weigths[k, ] <- ifelse(xtree$H1[k, ] > 0, AFw[i, k], 0)
                    }
                    i.primeW <- c(1:ncol(xtree.weigths))
                    for (k in 1:ncol(xtree.weigths)) {
                      if (sum(xtree.weigths[which(xtree.weigths[, k] > 
                                                  0), k]) != 0) {
                        i.primeW[k] <- mean(xtree.weigths[which(xtree.weigths[, 
                                                                              k] > 0), k], na.rm = TRUE)
                      } else {
                        i.primeW[k] <- 0
                      }
                    }
                    Out[i, 5] <- sum(i.primeW * xtree$h2.prime)
                    xtree.weigths <- xtree$H1
                    for (k in 1:nrow(S)) {
                      xtree.weigths[k, ] <- ifelse(xtree$H1[k, ] > 0, 
                                                   AFcomm[i, k], 0)
                    }
                    i.primeW <- c(1:ncol(xtree.weigths))
                    for (k in 1:ncol(xtree.weigths)) {
                      if (sum(xtree.weigths[which(xtree.weigths[, k] > 
                                                  0), k]) != 0) {
                        i.primeW[k] <- mean(xtree.weigths[which(xtree.weigths[, 
                                                                              k] > 0), k], na.rm = TRUE)
                      } else {
                        i.primeW[k] <- 0
                      }
                    }
                    Out[i, 6] <- sum(i.primeW * xtree$h2.prime)
                    return(Out)
                  }
  Out <- as.data.frame(na.omit(Out))
  Out <- Out[!duplicated(Out),]
  if (stand.FD == TRUE) {
    Out[, 4] <- Out[, 4]/max(Out[, 4])
  }
  row.names(Out) <- NULL
  Out
}
