###################################################################################################
# This tests the counting capabilities of various functions for DNase-C data.

Sys.setlocale(category="LC_COLLATE",locale="C")
dir.create("temp-dna")
file1 <- "temp-dna/1.h5"
file2 <- "temp-dna/2.h5"

suppressPackageStartupMessages(library(diffHic))

simDNA <- function(fout, chrs, npairs, rlen) { 
    r1 <- sample(length(chrs), npairs, replace=TRUE)
    r2 <- sample(length(chrs), npairs, replace=TRUE)
    p1 <- as.integer(runif(npairs, 1, chrs[r1] + 1))
    p2 <- as.integer(runif(npairs, 1, chrs[r2] + 1))
    l1 <- ifelse(rbinom(npairs, 1, 0.5)==1L, 1L, -1L)*rlen
    l2 <- ifelse(rbinom(npairs, 1, 0.5)==1L, 1L, -1L)*rlen

    savePairs(data.frame(anchor1.id=r1, anchor2.id=r2, anchor1.pos=p1, anchor2.pos=p2, anchor1.len=l1, anchor2.len=l2),
              fout, param=pairParam( GRanges(seqlengths=chrs) ))
    return(invisible(NULL))
}

comp <- function(chrs, npairs1, npairs2, dist, rlen=10, filter=1L, restrict=NULL, cap=NA) {
    simDNA(file1, chrs, npairs1, rlen)   
    simDNA(file2, chrs, npairs2, rlen)   

    # Output of squares.
    param <- pairParam( GRanges(seqlengths=chrs), restrict=restrict, cap=cap)
    y <- squareCounts(c(file1, file2), param=param, width=dist, filter=filter)

    # Reference. First, getting all bins.
    bin.coords <- list()
    for (i in names(chrs)) { 
        nbins <- ceiling(chrs[[i]]/dist)
        bin.ends <- pmin(chrs[[i]], seq_len(nbins)*dist)
        bin.starts <- c(1, head(bin.ends, -1)+1)
        bin.coords[[i]] <- GRanges(i, IRanges(bin.starts, bin.ends))
    }
    bin.offset <- c(0L, cumsum(lengths(bin.coords)))
    names(bin.coords) <- NULL
    suppressWarnings(bin.coords <- do.call(c, bin.coords))
    seqlengths(bin.coords) <- chrs
    bin.coords$nfrags <- 1L
    stopifnot(identical(bin.coords, regions(y)))

    # Now running through all bin pairs and assembling an InteractionSet.
    collected.isets <- list()
    collected.margins <- list()
    collected.totals <- list()
    for (f in 1:2) { 
        curf <- c(file1, file2)[f]
        fmat <- matrix(0L, length(bin.coords), length(bin.coords))
        total <- 0L

        for (i in seq_along(chrs)) { 
            for (j in seq_len(i)) {
                cur.i <- names(chrs)[i]
                cur.j <- names(chrs)[j]
                if (!is.null(restrict) && (!cur.i %in% restrict || !cur.j %in% restrict)) { next }
                curdat <- loadData(curf, cur.i, cur.j)
                total <- total+ nrow(curdat)
                
                p1 <- curdat$anchor1.pos + ifelse(curdat$anchor1.len > 0, 0L, -curdat$anchor1.len-1L)
                p1 <- pmin(p1, chrs[i])
                p2 <- curdat$anchor2.pos + ifelse(curdat$anchor2.len > 0, 0L, -curdat$anchor2.len-1L)
                p2 <- pmin(p2, chrs[j])
                b1 <- ceiling(p1/dist) + bin.offset[i]
                b2 <- ceiling(p2/dist) + bin.offset[j]

                for (x in seq_along(b1)) {
                    fmat[b1[x], b2[x]] <- fmat[b1[x], b2[x]] + 1L 
                    if (b1[x]!=b2[x]) { 
                        fmat[b2[x], b1[x]] <- fmat[b2[x], b1[x]] + 1L 
                    }
                }                
            }
        }

        cm <- ContactMatrix(fmat, seq_along(bin.coords), seq_along(bin.coords), regions=bin.coords)
        extractor <- upper.tri(fmat, diag=TRUE)
        is <- deflate(cm, extract=extractor)
        collected.isets[[f]] <- is
        collected.margins[[f]] <- as.integer(rowSums(fmat) + diag(fmat)) # diagonal gets counted twice.
        collected.totals[[f]] <- total
    }

    ref <- do.call(cbind, collected.isets)
    ref <- ref[rowSums(assay(ref)) >= filter,]
    interactions(ref) <- as(interactions(ref), "ReverseStrictGInteractions")
    storage.mode(assay(ref)) <- "integer"
    colnames(ref) <- NULL

    # Checking if interactions and counts are equal.
    m <- match(y, ref)
    stopifnot(identical(assay(ref)[m,], assay(y)))
    stopifnot(all(!is.na(m)))
    stopifnot(!anyDuplicated(m))
    stopifnot(nrow(y)==nrow(ref))
    
    # Checking the totals.
    stopifnot(identical(y$totals, as.integer(unlist(collected.totals))))
    if (is.null(restrict)) { 
        stopifnot(identical(y$totals, as.integer(c(npairs1, npairs2))))
    }
    totes <- totalCounts(c(file1, file2), param=param)
    stopifnot(identical(totes, y$totals))

    # Checking if the marginal counts... add up.
    mrg <- marginCounts(c(file1, file2), param=param, width=dist)
    out <- assay(mrg)
    dimnames(out) <- NULL
    stopifnot(identical(out, do.call(cbind, collected.margins)))
    stopifnot(identical(bin.coords, rowRanges(mrg)))

    # Checking that the neighborhood gives the same output.
    nbr <- neighborCounts(c(file1, file2), param=param, width=dist, filter=filter, flank=5)
    stopifnot(identical(assay(nbr), assay(y)))
    stopifnot(all(interactions(nbr)==interactions(y)))

    return(head(assay(y)))
}

chrs <- c(chrA=1000, chrB=2000)
comp(chrs, 100, 200, 100, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 200, 200, 100, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 100, 200, 75, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 100, 200, 220, rlen=10, filter=1L, restrict=NULL, cap=NA)

comp(chrs, 500, 500, 75, rlen=10, filter=5L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 100, rlen=10, filter=5L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 150, rlen=10, filter=5L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 200, rlen=10, filter=5L, restrict=NULL, cap=NA)

comp(chrs, 500, 500, 75, rlen=50, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 100, rlen=50, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 150, rlen=50, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 200, rlen=50, filter=1L, restrict=NULL, cap=NA)

comp(chrs, 500, 500, 75, rlen=10, filter=1L, restrict="chrA", cap=NA)
comp(chrs, 500, 500, 100, rlen=10, filter=1L, restrict="chrA", cap=NA)
comp(chrs, 500, 500, 150, rlen=10, filter=1L, restrict="chrA", cap=NA)
comp(chrs, 500, 500, 200, rlen=10, filter=1L, restrict="chrA", cap=NA)

comp(chrs, 500, 500, 75, rlen=10, filter=1L, restrict=NULL, cap=5) # Should have no effect.
comp(chrs, 500, 500, 100, rlen=10, filter=1L, restrict=NULL, cap=5)
comp(chrs, 500, 500, 150, rlen=10, filter=1L, restrict=NULL, cap=5)
comp(chrs, 500, 500, 200, rlen=10, filter=1L, restrict=NULL, cap=5)

chrs <- c(chrA=1000, chrB=1000, chrC=1000) # Trying for more chromosomes.
comp(chrs, 500, 500, 75, rlen=10, filter=1L, restrict=NULL, cap=NA) 
comp(chrs, 500, 500, 100, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 150, rlen=10, filter=1L, restrict=NULL, cap=NA)
comp(chrs, 500, 500, 200, rlen=10, filter=1L, restrict=NULL, cap=NA)

##################################################################################################
# Cleaning up.

unlink("temp-dna", recursive=TRUE)

##################################################################################################
# End.
