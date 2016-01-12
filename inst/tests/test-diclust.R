# This tests the post-hoc clustering methods. We're just recording results here,
# rather than doing rigorous tests, because that would be equivalent to just repeating
# the R code and that seems a bit like a waste of time. 

checkResults <- function(data.list, result.list, ..., true.pos) {
    out <- diClusters(data.list, result.list, ...)

    # Checking that the clustering is fine.
    all.ids <- unlist(out$id)
    ref <- do.call(c, data.list)
    bbox <- boundingBox(ref, all.ids)
    stopifnot(all(bbox$first==out$anchor1))
    stopifnot(all(bbox$second==out$anchor2))

    # Checking that the right interactions were chosen.
    all.ps <- sapply(result.list, FUN=function(x) { x$PValue })
    was.sig <- !is.na(all.ids)
    stopifnot(max(all.ps[was.sig]) > min(all.ps[!was.sig]))

    # Reporting the observed and estimated FDRs.
    np <- sum(!overlapsAny(GRangesList(out$anchor1, out$anchor2), true.pos))
    return(data.frame(Observed=np/length(out$anchor1), Estimated=out$FDR))
}

set.seed(100)
regions <- GRanges("chrA", IRanges(1:500, 1:500))
first.anchor <- sample(1000, 50, replace=TRUE)
second.anchor <- sample(1000, 50, replace=TRUE)
interactions <- InteractionSet(matrix(0, nrow=1000, ncol=1), ReverseStrictGInteractions(first.anchor, second.anchor, regions))
test.p <- runif(1000)
test.p[rep(1:2, 100) + rep(0:99, each=2) * 10] <- 0 

true.pos <- interactions[test.p==0]
checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=0, target=0.05, true.pos=true.pos)
checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=10, target=0.05, true.pos=true.pos)

checkResults(list(interactions, interactions[1:10]), list(data.frame(PValue=test.p), data.frame(PValue=test.p[1:10])), tol=0, target=0.05, true.pos=true.pos) # Multiple entries
checkResults(list(interactions, interactions[1:10]), list(data.frame(PValue=test.p), data.frame(PValue=test.p[1:10])), equiweight=FALSE, tol=0, target=0.05, true.pos=true.pos)

# Smaller number of DI entries
set.seed(50)
test.p <- runif(1000)
test.p[rep(1:2, 50) + rep(0:49, each=2) * 10] <- 0  

true.pos <- interactions[test.p==0]
checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=0, target=0.05, true.pos=true.pos)
checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=5, target=0.05, true.pos=true.pos)
checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=5, target=0.1, true.pos=true.pos)
checkResults(list(interactions), list(data.frame(whee=test.p)), tol=2, pval.col="whee", target=0.05, true.pos=true.pos)

signs <- c(-1, 1)[rbinom(100, 1, 0.5)+1]
checkResults(list(interactions, interactions[1:10]), list(data.frame(PValue=test.p, logFC=signs), data.frame(PValue=test.p[1:10], logFC=signs[1:10])), 
             tol=0, fc.col="logFC", target=0.05, true.pos=true.pos)

checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=0, grid.param=list(scale=5, iter=10), target=0.05, true.pos=true.pos) # Fiddling with grid search parameters.
checkResults(list(interactions), list(data.frame(PValue=test.p)), tol=0, grid.param=list(len=11, it=10), target=0.05, true.pos=true.pos)

###################################################
