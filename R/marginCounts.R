marginCounts <- function(files, param, width=50000)
# Gets the marginal counts i.e. sum of counts for each bin or region.
# This is useful to determine the `genomic coverage' of each region,
# based on the number of Hi-C read pairs involving that region.
#
# written by Aaron Lun
# some time ago.
# last modified 22 November 2015
{
	nlibs <- length(files)
	width <- as.integer(width)
	fragments <- param$fragments
	frag.by.chr <- .splitByChr(fragments)

	# Setting up other local references.
	restrict <- param$restrict
	discard <- .splitDiscards(param$discard)
	cap <- param$cap

	if (width < 0) { stop("width must be a non-negative integer") }
	new.pts <- .getBinID(fragments, width)
	bin.by.chr <- .splitByChr(new.pts$region)

	total.counts <- matrix(0L, length(new.pts$region), nlibs)
	full.sizes <- integer(nlibs)
	chrs <- seqlevels(fragments)

	# Running through each pair of chromosomes.
	overall <- .loadIndices(files, chrs, restrict)
	for (anchor in names(overall)) {
		current <- overall[[anchor]]
		first.anchor <- bin.by.chr$first[[anchor]]
		last.anchor <- bin.by.chr$last[[anchor]]
		nabins <- last.anchor - first.anchor + 1L
		keep.a <- first.anchor:last.anchor

		for (target in names(current)) {
			first.target <- bin.by.chr$first[[target]]
			last.target <- bin.by.chr$last[[target]]
			ntbins <- last.target - first.target + 1L
			keep.t <- first.target:last.target

			pairs <- .baseHiCParser(current[[target]], files, anchor, target,
				chr.limits=frag.by.chr, discard=discard, cap=cap)

			# Aggregating them for each library.
			for (lib in seq_len(nlibs)) {
				a.counts <- tabulate(new.pts$id[pairs[[lib]]$anchor.id]-first.anchor+1L, nbins=nabins)
				total.counts[keep.a,lib] <- total.counts[keep.a,lib] + a.counts
				t.counts <- tabulate(new.pts$id[pairs[[lib]]$target.id]-first.target+1L, nbins=ntbins)
				total.counts[keep.t,lib] <- total.counts[keep.t,lib] + t.counts
				full.sizes[lib] <- full.sizes[lib] + nrow(pairs[[lib]])
			}
		}	
	}
	
	# Aggregating all elements.
	return(RangedSummarizedExperiment(list(counts=total.counts), colData=DataFrame(totals=full.sizes), 
			rowRanges=new.pts$region, metadata=List(param=param)))
}

