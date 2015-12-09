prepPseudoPairs <- function(bam, param, file, dedup=TRUE, ichim=TRUE, chim.span=1000, minq=NA, output.dir=NULL)
# This function acts the same as preparePairs, but it assumes that you're
# putting things into contiguous bins across the genome. The idea is to
# allow DNase-digested Hi-C experiments to fit in the pipeline, where reads
# are assigned to pseudo-restriction fragments.
#
# written by Aaron Lun
# created 27 March 2015
# last modified 9 December 2015
{
	fragments <- param$fragments
	n.per.chr <- runLength(seqnames(fragments))
	frag.data <- .splitByChr(fragments)

	chrs <- frag.data$chr
	last.in.chr <- frag.data$last
	before.first <- as.list(c(0L, last.in.chr[-length(chrs)]))
	names(before.first) <- chrs

	bin.width <- max(width(fragments))
	if (!all(bin.width==width(fragments)[-last.in.chr])) {
		stop("pseudo-fragments should be constant size") 
	}

	# Checking consistency between SAM chromosome lengths and the ones in the cuts.
	chromosomes<-scanBamHeader(bam)[[1]]$targets
    ref.chrs <- names(chromosomes)
    m <- match(ref.chrs, chrs)
    if (any(is.na(m))) { stop("missing chromosomes in cut site list") }

    n.per.chr <- n.per.chr[m] # Rearranging so that BAM alignment TIDs refer to the correct chromosomes.
    endpoints <- end(fragments)[last.in.chr[m]]   
	for (x in seq_along(ref.chrs)) {
		if (chromosomes[x]!=endpoints[x]) {
			stop("length of ", ref.chrs[x], " is not consistent between BAM file and fragment ranges")
		}
	}

	# Enforcing input types.
	minq <- as.integer(minq)
	ichim <- as.logical(ichim)
	chim.span <- as.integer(chim.span)
	dedup <- as.logical(dedup)

    # Setting up the output directory.
    if (is.null(output.dir)) { 
        output.dir <- tempfile(tmpdir=".")
    }
    if (file.exists(output.dir)) { 
        stop("output directory already exists")
    }
    dir.create(output.dir)
    on.exit({ unlink(output.dir, recursive=TRUE) })
    prefix <- file.path(output.dir, "")

    out <- .Call(cxx_report_hic_binned_pairs, n.per.chr, bin.width, bam, prefix, !ichim, chim.span, minq, dedup)
    if (is.character(out)) { stop(out) }
    .process_output(out, file, ref.chrs, before.first)
}

segmentGenome <- function(bs, size=500) 
# Segments the genome into pseudo-fragments, for use in assigning reads 
# from a DNase Hi-C experiment.
#
# written by Aaron Lun
# created 27 March 2015
{
	if (is(bs, "BSgenome")) {
		ref.len <- seqlengths(bs)
		gen <- genome(bs)
	} else if (is.character(bs)) {
		bs <- readDNAStringSet(bs)
		ref.len <- width(bs)
		names(ref.len) <- names(bs)
		gen <- NA
	} else {
		ref.len <- bs
		gen <- NA
	}

	everything <- list()
	for (chr in names(ref.len)) {
		bin.dex <- seq_len(ceiling(ref.len[[chr]]/size))
		current <- GRanges(chr, IRanges((bin.dex - 1L)*size + 1L, bin.dex*size))
		if (end(current)[length(current)] > ref.len[[chr]]) { 
			end(current)[length(current)] <- ref.len[[chr]]
		}
		everything[[length(everything)+1L]] <- current
	}
	suppressWarnings(everything <- do.call(c, everything))
	seqlengths(everything) <- ref.len
	return(everything)
}
