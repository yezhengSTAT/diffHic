.splitByChr <- function(ranges)
# Gets the start and end for each chromosome in the sorted GRanges. 
{
    chrs <- as.character(runValue(seqnames(ranges)))
    if (anyDuplicated(chrs)) { stop("ranges for each chromosome should be consecutive") }
    ref.len <- runLength(seqnames(ranges))
    end.index <- cumsum(ref.len)
    start.index <- end.index - ref.len + 1L
    names(end.index) <- names(start.index) <- chrs
    return(list(chr=chrs, first=start.index, last=end.index))
}

####################################################################################################

.process_output <- function(c_out, file, chrs, chr.start) 
# Converts the output of the C++ code in preparePairs or prepPseudoPairs
# into HDF5 files. Also formats and returns the diagnostics.
{
    .initializeH5(file)
    for (a1.dex in seq_along(c_out[[1]])) { 
        curnames <- c_out[[1]][[a1.dex]]
        not.empty <- curnames!=""
        if (!any(not.empty)) { next }
        anchor1 <- chrs[a1.dex]
        .addGroup(file, anchor1)

        for (a2.dex in which(not.empty)) { 
            anchor2 <- chrs[a2.dex]
            current.file <- curnames[a2.dex]
            out <- read.table(current.file, header=FALSE, colClasses="integer")
            colnames(out) <- c("anchor1.id", "anchor2.id", "anchor1.pos", "anchor2.pos", "anchor1.len", "anchor2.len")

            out$anchor1.id <- out$anchor1.id+chr.start[[anchor1]]
            out$anchor2.id <- out$anchor2.id+chr.start[[anchor2]]
            out <- out[order(out$anchor1.id, out$anchor2.id),,drop=FALSE]
            rownames(out)<-NULL
            .writePairs(out, file, anchor1, anchor2)
        }
    }

    c_out <- c_out[-1]
    names(c_out) <- c("pairs", "same.id", "singles", "chimeras")
    names(c_out$pairs) <-c("total", "marked", "filtered", "mapped")
    names(c_out$same.id) <- c("dangling", "self.circle")
    names(c_out$chimeras) <- c("total", "mapped", "multi", "invalid")
    return(c_out)
}

####################################################################################################

.getBinID <- function(fragments, width) 
# Determines which bin each restriction fragment is in. Also records the rounded
# start and stop site for each bin. Returns a set of bin ids for each restriction
# fragment on each chromosome, as well as the coordinates of each bin.
{
    width<-as.integer(width)
    if (length(fragments) == 0L) {
        return(.createBins(fragments, width))
    }

    out.ids<-integer(length(fragments))
    out.ranges<-list()
    last<-0L
    frag.data <- .splitByChr(fragments)
    nfrags <- list() 
    
    for (x in seq_along(frag.data$chr)) {
        curindex <- frag.data$first[x]:frag.data$last[x]
        curf <- fragments[curindex]
        mids <- (start(curf)+end(curf))/2
        bin.id <- as.integer((mids-0.1)/width)+1L 
        # The '-0.1' in the preceding step reduces 'mids' that are exact multiples 
        # of 'width', so each bin is from (n*width, (n+1)*width] for integer 'n'.

        processed <- rle(bin.id)
        ns <- length(processed$value)
        processed$values <- seq_len(ns)
        nfrags[[x]] <- processed$length
        out.ids[curindex] <- inverse.rle(processed)+last
        
        endx <- cumsum(processed$length)
        startx <- rep(1L, ns)
        if (ns>=2L) { startx[-1] <- endx[-ns]+1L }
        out.ranges[[x]] <- GRanges(frag.data$chr[x], IRanges(start(curf[startx]), end(curf[endx])))
        last <- last+ns
    }

    # Wrapping up.
    suppressWarnings(out.ranges <- do.call(c, out.ranges))
    seqlevels(out.ranges) <- seqlevels(fragments)
    seqlengths(out.ranges) <- seqlengths(fragments)
    out.ranges$nfrags <- unlist(nfrags)
    return(list(id=out.ids, region=out.ranges))
}

.createBins <- function(fragments, width) 
# This creates regular contiguous bins of size 'width'. Each bin
# is assigned to itself; allocation of read pairs into bins is done below.
# This allows free-floating bins for use with DNase-C data.
{
    ref.len <- seqlengths(fragments)
    everything <- list()
    for (chr in names(ref.len)) {
        bin.dex <- seq_len(ceiling(ref.len[[chr]]/width))
        end.pt <- pmin(bin.dex * width, ref.len[[chr]])
        current <- GRanges(chr, IRanges((bin.dex - 1L)*width + 1L, end.pt))
        everything[[length(everything)+1L]] <- current
    }
    suppressWarnings(everything <- do.call(c, everything))
    seqlengths(everything) <- ref.len
    return(list(id=seq_along(everything), region=everything))
}

####################################################################################################

.parseParam <- function(param) {
    output <- list()
    
    fragments <- param$fragments
    if (length(fragments)==0L) { 
        output$chrs <- seqlevels(fragments)
        output$frag.by.chr <- list(chr=chrs, 
                                   first=setNames(integer(
 
    } else {
        output$chrs <- seqlevelsInUse(fragments)
        output$frag.by.chr <- .splitByChr(fragments)
    }
}

####################################################################################################

.baseHiCParser <- function(ok, files, anchor1, anchor2, chr.limits, discard, cap, width=NA)
# A convenience function for loading counts from file for a given anchor/anchor pair.
# It will also bin the read pairs if 'width' is specified (for DNase-C experiments).
{
    overall<-list()
    adisc <- discard[[anchor1]]
    tdisc <- discard[[anchor2]]
    do.cap <- !is.na(cap) 
    do.bin <- !is.na(width)

    for (x in seq_along(ok)) {
        if (!ok[x]) { 
            overall[[x]] <- data.frame(anchor1.id=integer(0), anchor2.id=integer(0))
        } else {
            out <- .getPairs(files[x], anchor1, anchor2)
    
            # Checking fidelity of the input.
            check <- .Call(cxx_check_input, out$anchor1.id, out$anchor2.id)
            if (is.character(check)) { stop(check) }

            # Checking that we're all on the right chromosome.
            if (nrow(out)) { 
                if (max(out$anchor1.id) > chr.limits$last[[anchor1]] || 
                        min(out$anchor1.id) < chr.limits$first[[anchor1]]) { 
                    stop("anchor1 index outside range of fragment object") 
                }
                if (max(out$anchor2.id) > chr.limits$last[[anchor2]] || 
                        min(out$anchor2.id) < chr.limits$first[[anchor2]]) { 
                    stop("anchor2 index outside range of fragment object") 
                }
            }

            # Overlapping with those in the discard intervals.
            if (!is.null(adisc) || !is.null(tdisc)) {
                a.hits <- t.hits <- FALSE
                if (!is.null(adisc)) {
                    a.hits <- overlapsAny(IRanges(out$anchor1.pos, out$anchor1.pos+abs(out$anchor1.len)-1L), adisc, type="within")
                }
                if (!is.null(tdisc)) { 
                    t.hits <- overlapsAny(IRanges(out$anchor2.pos, out$anchor2.pos+abs(out$anchor2.len)-1L), tdisc, type="within")
                }
                out <- out[!a.hits & !t.hits,,drop=FALSE]
            }

            # Removing read pairs above the cap for each restriction fragment pair.
            if (do.cap) { 
                capped <- .Call(cxx_cap_input, out$anchor1.id, out$anchor2.id, cap)
                if (is.character(capped)) { stop(capped) }
                out <- out[capped,]
            }

            if (!is.na(width)) { out <- .binReads(out, width) } 
            dim(out$anchor1.id) <- dim(out$anchor2.id) <- NULL
            overall[[x]] <- out[,c("anchor1.id", "anchor2.id")]
        }
    }
    return(overall)
}

.binReads <- function(pairs, width) {
    a1.5pos <- pairs$anchor1.pos
    a1.5len <- pairs$anchor1.len
    a1.r <- a1.5len < 0L
    a1.5pos[a1.r] <- a1.5pos[a1.r] + a1.len[a1.r] - 1L

    a2.5pos <- pairs$anchor2.pos
    a2.5len <- pairs$anchor2.len
    a2.r <- a2.5len < 0L
    a2.5pos[a2.r] <- a2.5pos[a2.r] + a2.len[a2.r] - 1L

    pairs$anchor1.id <- ceiling(a1.5pos/width)
    pairs$anchor2.id <- ceiling(a2.5pos/width)
    o <- order(pairs$anchor1.id, pairs$anchor2.id)
    pairs[o,]
}

####################################################################################################

.splitDiscards <- function(discard) 
# Splits the discard GRanges into a list of constituent chromosomes,
# along with IRanges for everything. This allows easy access to the
# ranges on individual chromosomes, rather than overlapping with everything.
{
    if (is.null(discard) || length(discard)==0L) { return(NULL) }
    discard <- sort(discard)
    all.chrs <- as.character(runValue(seqnames(discard)))
    all.len <- runLength(seqnames(discard))
    chr.ends <- cumsum(all.len)
    chr.starts <- c(1L, chr.ends[-length(chr.ends)]+1L)

    output <- list()
    for (i in seq_along(all.chrs)) {
        chr <- all.chrs[i]
        ix <- chr.starts[i]:chr.ends[i]
        output[[chr]] <- reduce(ranges(discard[ix]))
    }

    return(output)
}

####################################################################################################

.getChrsInUse <- function(fragments) 
# Gets the chromosomes that are available. If we're 
# working with DNase-C data, this is all seqlevels,
# otherwise it's just the ones that are in use.    
{
}
