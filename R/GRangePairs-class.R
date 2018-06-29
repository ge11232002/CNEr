# This GRangePairs directly extends `Pairs` class.
### -----------------------------------------------------------------
### GRangePairs: class
### Exported!
setClass(Class="GRangePairs",
         contains="Pairs",
         slots=c(first="GRanges",
                 second="GRanges"),
         prototype=list(first=GRanges(),
                        second=GRanges()
                        ))

setValidity("GRangePairs",
            function(object){
              x_first <- object@first
              x_second <- object@second
              ## Test first's class
              if(class(x_first) != "GRanges")
                return("'x@first' must be a GRanges instance")
              ## test second's class
              if(class(x_second) != "GRanges")
                return("'x@second' must be a GRanges instance")
              return(TRUE)
            })

### Formal API:
###   GRangePairs(x) - constructor.
###   names(x)    - NULL or character vector.
###   length(x)   - single integer N. Nb of pairs in 'x'.
###   first(x)    - returns "first" slot.
###   last(x)     - returns "last" slot.
###   seqnames(x) - returns DataFrame of seqnames of first, last GRanges.
###   strand(x)   - returns DataFrame of strands of first, last GRanges.
###   seqinfo(x)  - returns list of seqinfo of first, last GRanges.
###   x[i]        - GRangePairs object of the same class as 'x'
###   x[[i]]      - GRanges object of concatenating the i-th GRangesPairs's
###                   first, last GRanges.
###   unlist(x)   - unlist the x into a GRanges object by concatenating
###                   each pair first.
###   grglist(x)  - GRangesList object of the same length as 'x'.
###   show(x)     - compact display in a data.frame-like fashion.

### -----------------------------------------------------------------
### GRangePairs Constructor.
### Exported!
GRangePairs <- function(first=GRanges(), second=GRanges(), ..., names=NULL,
                        hits=NULL){
  if(!is.null(hits)) {
    stopifnot(is(hits, "Hits"),
              queryLength(hits) == length(first),
              subjectLength(hits) == length(second))
    first <- first[queryHits(hits)]
    second <- second[subjectHits(hits)]
  }
  if(!(is(first, "GRanges") && is(second, "GRanges")))
    stop("'first' and 'second' must be GRanges objects")
  stopifnot(length(first) == length(second),
            is.null(names) || length(names) == length(first))
  if(!missing(...)) {
    elementMetadata <- DataFrame(...)
  }else{
    elementMetadata <- S4Vectors:::make_zero_col_DataFrame(length(first))
  }
  #rownames(elementMetadata) <- names
  new("GRangePairs",first=first, second=second, NAMES=names,
      elementMetadata=elementMetadata)
}

### -----------------------------------------------------------------
### GRangePairs getters and setters
### Exported!
setMethod("last", "GRangePairs",
          function(x)
          {
            second(x)
          }
)

setMethod("seqnames", "GRangePairs",
          function(x){
            ans <- DataFrame(first=seqnames(first(x)),
                             second=seqnames(second(x)))
            ans
          }
)

setMethod("strand", "GRangePairs",
          function(x){
            ans <- DataFrame(first=strand(first(x)),
                             second=strand(second(x)))
            ans
          }
)

setMethod("seqinfo", "GRangePairs",
          function(x) list(seqinfoFirst=seqinfo(first(x)),
                           seqinfoSecond=seqinfo(second(x)))
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combine
### Exported!
.unlist_list_of_GRangePairs <- function(x, Class){
  metadata <- do.call(rbind, lapply(x, mcols))
  rownames(metadata) <- NULL
  new(Class, first=do.call(c, lapply(x, first)),
      second=do.call(c, lapply(x, second)),
      elementMetadata=metadata,
      ### FIXME: breaks if only some names are NULL
      NAMES=unlist(lapply(x, names)))
}
setMethod("c", "GRangePairs",
          function(x, ..., recursive=FALSE){
            if(isTRUE(recursive))
              stop("'recursive' argument not supported")
            if (missing(x))
              args <- unname(list(...))
            else args <- unname(list(x, ...))
            .unlist_list_of_GRangePairs(args, class(args[[1]]))
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List methods.
###

### TODO: Remove this method after the definition of the GAlignmentPairs
### class is changed to derive from CompressedList.
setMethod("unlist", "GRangePairs",
          function(x, use.names=TRUE)
          {
            if (!isTRUEorFALSE(use.names))
              stop("'use.names' must be TRUE or FALSE")
            x_first <- first(x)
            x_last <- second(x)
            collate_subscript <-
              S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(x))
            ans <- c(x_first, x_last)[collate_subscript]
            if (use.names)
              names(ans) <- rep(names(x), each=2L)
            ans
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Coercion.
###
setMethod("grglist", "GRangePairs",
          function(x, use.mcols=FALSE)
          {
            if (!isTRUEorFALSE(use.mcols))
              stop("'use.mcols' must be TRUE or FALSE")
            x_mcols <- mcols(x)
            if (use.mcols && "query.break" %in% colnames(x_mcols))
              stop("'mcols(x)' cannot have reserved column \"query.break\"")
            x_first <- first(x)
            x_last <- second(x)
            collate_subscript <-
              S4Vectors:::make_XYZxyz_to_XxYyZz_subscript(length(x))
            x_unlisted <- c(x_first, x_last)
            x_unlisted <- x_unlisted[collate_subscript]
            grl <- as(x_unlisted, "GRangesList")
            ans <- GenomicAlignments:::shrinkByHalf(grl)
            names(ans) <- names(x)
            ans_mcols <- DataFrame(query.break=mcols(ans)$nelt1)
            if (use.mcols)
              ans_mcols <- cbind(ans_mcols, x_mcols)
            mcols(ans) <- ans_mcols
            ans
          }
)

setAs("GRangePairs", "GRangesList",
      function(from) grglist(from, use.mcols=TRUE)
)
setAs("GRangePairs", "GRanges",
      function(from) unlist(from, use.names=TRUE)
)

### -----------------------------------------------------------------
### swap method for GRangePairs: first becomes last and last becomes first
### Exported!
setGeneric("swap", function(x) standardGeneric("swap"))
setMethod("swap", "GRangePairs", function(x){
  BiocGenerics:::replaceSlots(x, first=second(x), second=first(x))
})

### -----------------------------------------------------------------
### unique: keep the unique GRangePairs
### Exported!
setMethod("unique", "GRangePairs", function(x){
  duplicatedIndex <- duplicated(paste(paste(first(x)), paste(second(x))))
  x[!duplicatedIndex]
})