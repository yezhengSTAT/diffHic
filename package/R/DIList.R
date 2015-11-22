# Defines the DIList class, that will be used to hold various information.

setClass("DIList", representation(counts="matrix", 
			anchors="integer", targets="integer", regions="GRanges",
			colData="DataFrame", exptData="List"))

setValidity("DIList", function(object) {
	if (nrow(object@counts)!=length(object@anchors)) {
		return('rows in count matrix not equal to length of anchor vector')
	} 
	if (nrow(object@counts)!=length(object@targets)) { 
		return('rows in count matrix not equal to length of target vector')
	}
	if (ncol(object@counts)!=nrow(object@colData)) {
		return('columns of count matrix not equal to rows of column data frame')
	}

	if (!all(object@anchors >= 1L)) { 
		return('not all anchors are positive integers')
	} 
	if (!all(object@targets >= 1L)) {
		return('not all targets are positive integers')
	}
	if (!all(object@anchors <= length(object@regions))) {
		return('not all anchors refer to valid regions')
	} 
	if (!all(object@targets <= length(object@regions))) { 
		return('not all targets refer to valid regions')
	}
	if (!all(object@anchors >= object@targets)) { 
		return('target indices cannot be greater than anchor indices')
	}
	return(TRUE)
})

setMethod("initialize", signature("DIList"), function(.Object, ...) {
	value <- callNextMethod()
	validObject(value)
	value
})

setMethod("show", signature("DIList"), function(object) {
	total <- nrow(object@counts)
	nregs <- length(object@regions)
	nlibs <- ncol(object@counts)
	cat("DIList object for", nlibs, ifelse(nlibs==1L, "library", "libraries"), 
		"with", total, ifelse(total==1L, "pair", "pairs"), "across", 
		nregs, ifelse(nregs==1L, "region\n", "regions\n"))
})

DI2IS <- function(x) {
    InteractionSet(list(counts=x@counts), anchor1=x@anchor1, anchor2=x@anchor2,
                   colData=x@colData, regions=x@regions, metadata=x@exptData)
}
