# Tests for zero-inputs into various diffHic functions.

suppressWarnings(suppressPackageStartupMessages(require(diffHic)))

# Don't worry about cases involving empty files; 
# that just shouldn't be a practical concern.

ghost <- InteractionSet(matrix(0, nrow=0, ncol=1), 
    GInteractions(integer(0), integer(0), regions=GRanges("chrA", IRanges(1:5, 1:5))),
    colData=DataFrame(totals=0))

getArea(ghost)
getArea(ghost, bp=TRUE)

f.out <- "empty.h5"
param <- pairParam(GRanges("chrA", IRanges(1:5, 1:5)))
savePairs(data.frame(anchor1.id=integer(0), anchor2.id=integer(0)), file=f.out, param=param)

loadChromos(f.out) # While we're here, let's see what happens.
unlink(f.out)

filterDirect(ghost)
filterDirect(ghost, reference=ghost)

filterTrended(ghost)
try(filterTrended(ghost, reference=ghost)) # This will fail, as interpolation is impossible.

filterPeaks(ghost, integer(0))

enrichedPairs(ghost, abundances=numeric(0))

try(compartmentalize(ghost)) # This will fail, as interpolation is impossible.
try(compartmentalize(ghost, dist.correct=FALSE)) # This will also fail, due to non-unique k-means.

correctedContact(ghost)

ghost.ranges <- SummarizedExperiment(matrix(0, 0, 0), GRanges())
try(normalizeCNV(ghost, ghost.ranges)) # locfit isn't as robust as loessFit

matchMargins(ghost, ghost.ranges)

asDGEList(ghost)

normOffsets(ghost)
normOffsets(ghost, type="loess")

