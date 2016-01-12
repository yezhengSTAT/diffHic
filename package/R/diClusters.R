diClusters <- function(data.list, result.list, target, equiweight=TRUE, cluster.args=list(), pval.col="PValue", fc.col=NA, grid.param=NULL) 
# Performs post-hoc clustering of significant bin pairs, to control the 
# cluster-level FDR at 'target'. Adjusts the weights so that the contribution
# from each set of bin pairs is the same.
#
# written by Aaron Lun
# created 12 January 2016
{
    # Setting initial parameters.    
    if (missing(target)) {
        target <- 0.05
        warning("unspecified 'target' for the cluster-level FDR set to 0.05")
    }
    expanded.args <- as.list(match.call(clusterPairs, do.call(call, c("clusterPairs", cluster.args))))
    if (is.null(expanded.args$tol)) {
        cluster.args$tol <- 1
        warning("'tol' for 'clusterPairs' set to a default of 1 bp")
    }

    # Checking inputs.
	nset <- length(data.list)
	if (nset!=length(result.list)) { stop("indices must have same length as result list") }
	for (x in seq_len(nset)) {
		if (!identical(length(indices[[x]]), nrow(result.list[[x]]))) {
 		   	stop("corresponding entries of data and result lists must have same number of entries") 
        }
	}

    # Computing the adjusted FDR for all of these samples (with equi-weighting).
    if (equiweight) {
        weights <- rep(1/lengths(data.list), lengths(data.list))
    } else {
        weights <- NULL
    }

    all.ps <- list()
    for (x in seq_len(nset)) { 
        all.ps[[x]] <- result.list[[x]][,pval.col]
    }
    all.ps <- unlist(all.ps)
    adjp <- csaw:::.weightedFDR(all.ps, weights)

    in.each.group <- list()
    last <- 0L
    for (x in seq_len(nset)) { 
        in.each.group[[x]] <- last + seq_len(nrow(result.list[[x]]))
        last <- last + nrow(result.list[[x]])
    }

    # Getting the sign.
    all.signs <- list()
    if (is.na(fc.col)) { 
        all.signs <- lapply(lengths(data.list), logical)
    } else {
        all.signs <- lapply(data.list, function(x) { x[,fc.col] > 0 })
    }

	# Controlling the cluster-level FDR.
    FUN <- function(sig) {
        pos.data.list <- neg.data.list <- list()
        for (x in seq_len(nset)) { 
            cur.sig <- sig[in.each.group[[x]]]
            pos.data.list[[x]] <- data.list[[x]][cur.sig & all.signs[[x]],]
            neg.data.list[[x]] <- data.list[[x]][cur.sig & !all.signs[[x]],]
        }
        pos.clust <- do.call(clusterPairs, c(neg.data.list, cluster.args))
        neg.clust <- do.call(clusterPairs, c(pos.data.list, cluster.args))
    
        # Assembling it back into a single return value.
        clust.indices <- list()
        additional <- length(pos.clust$anchor1)
        for (x in seq_len(nset)) {
            cur.sig <- sig[in.each.group[[x]]]
            cur.signs <- all.signs[[x]][cur.sig] 
            full.ids <- as.integer(cur.sig)
            full.ids[cur.signs] <- pos.clust$indices[[x]]
            full.ids[!cur.signs] <- neg.clust$indices[[x]] + additional
            clust.indices[[x]] <- full.ids
        }
        names(clust.indices) <- names(data.list)
        list(indices=clust.indices, anchor1=c(pos.clust$anchor1, neg.clust$anchor1),
                anchor2=c(pos.clust$anchor2, neg.clust$anchor2))
    }
    out <- controlClusterFDR(target=target, adjp=adjp, FUN=function(sig) { unlist(FUN(sig)$indices) },
                             weight=weights, grid.param=grid.param)
    sig <- adjp <= out$threshold
    clusters <- FUN(sig)

    # Cleaning up the output.
    sig.by.group <- split(sig, groupings)
    for (x in seq_len(nset)) {
        full.ids <- rep(NA_integer_, data.list[[x]])
        full.ids[sig.by.group[[x]]] <- clusters$indices[[x]]
        clusters$indices[[x]] <- full.ids
    }
    clusters$FDR <- out$FDR    
    return(clusters)
}
