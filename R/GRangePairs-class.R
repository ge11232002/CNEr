# This GRangePairs is purely based on the implementation of GAlignmentPairs
### -----------------------------------------------------------------
### GRangePairs: class
### Exported!
setClass(Class="GRangePairs",
         contains="List",
         slots=c(NAMES="characterORNULL",      # R doesn't like @names !!
                 first="GRanges",
                 last="GRanges",
                 elementMetadata="DataFrame"),
         prototype=list(elementType="GRanges",
                        elementMetadata=DataFrame()))

setValidity("GRangePairs",
            function(object){
              x_first <- object@first
              x_last <- object@last
              ## Test NAMES
              x_names <- object@NAMES
              if (is.null(x_names))
                return(NULL)
              if (!is.character(x_names) || !is.null(attributes(x_names))) {
                msg <- c("'names(x)' must be NULL or a character vector ",
                         "with no attributes")
                return(paste(msg, collapse=""))
              }
              if (length(x_names) != length(object))
                return("'names(x)' and 'x' must have the same length")
              ## Test first's class
              if(class(x_first) != "GRanges")
                return("'x@first' must be a GRanges instance")
              ## test last's class
              if(class(x_last) != "GRanges")
                return("'x@last' must be a GRanges instance")
              ## test the size of first and last
              if (length(x_last) != length(x_first))
                return("'x@last' and 'x@first' must have the same length")
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
GRangePairs <- function(first=GRanges(), last=GRanges(), names=NULL){
  if (!(is(first, "GRanges") && is(last, "GRanges")))
    stop("'first' and 'last' must be GRanges objects")
  if (length(first) != length(last))
    stop("'first' and 'last' must have the same length")
  new("GRangePairs", NAMES=names, first=first, last=last,
      elementMetadata=S4Vectors:::make_zero_col_DataFrame(length(first)))
}


### -----------------------------------------------------------------
### GRangePairs getters and setters
### Exported!
setMethod("names", "GRangePairs", function(x) x@NAMES)
setMethod("length", "GRangePairs", function(x) length(x@first))
setMethod("first", "GRangePairs",
          function(x, real.strand=FALSE, invert.strand=FALSE)
          {
            ans <- setNames(x@first, names(x))
            ans
          }
)

setMethod("last", "GRangePairs",
          function(x, real.strand=FALSE, invert.strand=FALSE)
          {
            ans <- setNames(x@last, names(x))
            ans
          }
)

setMethod("seqnames", "GRangePairs",
          function(x){
            ans <- DataFrame(first=seqnames(x@first),
                             last=seqnames(x@last))
            ans
          }
)

setMethod("strand", "GRangePairs",
          function(x){
            ans <- DataFrame(first=strand(x@first),
                             last=strand(x@last))
            ans
          }
)

setMethod("seqinfo", "GRangePairs",
          function(x) list(seqinfoFirst=seqinfo(x@first),
                           seqinfoLast=seqinfo(x@last))
)

setReplaceMethod("names", "GRangePairs",
                 function(x, value)
                 {
                   if (!is.null(value))
                     value <- as.character(value)
                   x@NAMES <- value
                   validObject(x)
                   x
                 }
)

### -----------------------------------------------------------------
### Vector methods.
###
setMethod("extractROWS", "GRangePairs",
          function(x, i)
          {
            i <- normalizeSingleBracketSubscript(i, x, as.NSBS=TRUE)
            ans_NAMES <- extractROWS(x@NAMES, i)
            ans_first <- extractROWS(x@first, i)
            ans_last <- extractROWS(x@last, i)
            ans_elementMetadata <- extractROWS(x@elementMetadata, i)
            BiocGenerics:::replaceSlots(x,
                                        NAMES=ans_NAMES,
                                        first=ans_first,
                                        last=ans_last,
                                        elementMetadata=ans_elementMetadata)
          }
)

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### List methods.
###
.GRangePairs.getElement <- function(x, i)
{
  c(x@first[i], x@last[i])
}

setMethod("[[", "GRangePairs",
          function(x, i, j, ... , drop=TRUE)
          {
            if (missing(i) || !missing(j) || length(list(...)) > 0L)
              stop("invalid subsetting")
            i <- normalizeDoubleBracketSubscript(i, x)
            .GRangePairs.getElement(x, i)
          }
)

### TODO: Remove this method after the definition of the GAlignmentPairs
### class is changed to derive from CompressedList.
setMethod("unlist", "GRangePairs",
          function(x, use.names=TRUE)
          {
            if (!isTRUEorFALSE(use.names))
              stop("'use.names' must be TRUE or FALSE")
            x_first <- x@first
            x_last <- x@last
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
            x_first <- x@first
            x_last <- x@last
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
setAs("GRangePairs", "DataFrame", function(from) {
  firstDF <- as(first(from), "DataFrame")
  colnames(firstDF) <- paste0(colnames(firstDF), ".first")
  lastDF <- as(last(from), "DataFrame")
  colnames(lastDF) <- paste0(colnames(lastDF), ".last")
  DataFrame(firstDF, lastDF, mcols(from))
})

setMethod("as.data.frame", "GRangePairs",
          function(x, row.names = NULL, optional = FALSE) {
            as.data.frame(as(x, "DataFrame"), row.names=row.names,
                          optional=optional)
          })

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Combining.
###
### 'Class' must be "GRangePairs" or the name of a concrete subclass of
### GRangePairs
### 'objects' must be a list of GRangePairs objects.
### Returns an instance of class 'Class'.
combine_GRangePairs_objects <- function(Class, objects,
                                        use.names=TRUE, ignore.mcols=FALSE)
{
  if (!isSingleString(Class))
    stop("'Class' must be a single character string")
  if (!extends(Class, "GRangePairs"))
    stop("'Class' must be the name of a class that extends GRangePairs")
  if (!is.list(objects))
    stop("'objects' must be a list")
  if (!isTRUEorFALSE(use.names))
    stop("'use.names' must be TRUE or FALSE")
  ### TODO: Support 'use.names=TRUE'.
  if (use.names)
    stop("'use.names=TRUE' is not supported yet")
  if (!isTRUEorFALSE(ignore.mcols))
    stop("'ignore.mcols' must be TRUE or FALSE")

  if (length(objects) != 0L) {
    ## TODO: Implement (in C) fast 'elementIsNull(objects)' in IRanges,
    ## that does 'sapply(objects, is.null, USE.NAMES=FALSE)', and use it
    ## here.
    null_idx <- which(sapply(objects, is.null, USE.NAMES=FALSE))
    if (length(null_idx) != 0L)
      objects <- objects[-null_idx]
  }
  if (length(objects) == 0L)
    return(new(Class))
  ## TODO: Implement (in C) fast 'elementIs(objects, class)' in IRanges, that
  ## does 'sapply(objects, is, class, USE.NAMES=FALSE)', and use it here.
  ## 'elementIs(objects, "NULL")' should work and be equivalent to
  ## 'elementIsNull(objects)'.
  if (!all(sapply(objects, is, Class, USE.NAMES=FALSE)))
    stop("the objects to combine must be ", Class, " objects (or NULLs)")
  objects_names <- names(objects)
  names(objects) <- NULL  # so lapply(objects, ...) below returns an
  # unnamed list

  ## Combine "NAMES" slots.
  NAMES_slots <- lapply(objects, function(x) x@NAMES)
  ## TODO: Use elementIsNull() here when it becomes available.
  has_no_names <- sapply(NAMES_slots, is.null, USE.NAMES=FALSE)
  if (all(has_no_names)) {
    ans_NAMES <- NULL
  } else {
    noname_idx <- which(has_no_names)
    if (length(noname_idx) != 0L)
      NAMES_slots[noname_idx] <-
        lapply(elementNROWS(objects[noname_idx]), character)
    ans_NAMES <- unlist(NAMES_slots, use.names=FALSE)
  }

  ## Combine "first" slots.
  first_slots <- lapply(objects, function(x) x@first)
  ans_first <- do.call(c, c(first_slots, ignore.mcols=ignore.mcols))
  
  ## Combine "last" slots.
  last_slots <- lapply(objects, function(x) x@last)
  ans_last <- do.call(c, c(last_slots, ignore.mcols=ignore.mcols))
  
  ## Combine "mcols" slots. We don't need to use fancy
  ## IRanges:::rbind.mcols() for this because the "mcols" slot of a
  ## GAlignmentPairs object is guaranteed to be a DataFrame.
  if (ignore.mcols) {
    ans_mcols <- S4Vectors:::make_zero_col_DataFrame(length(ans_first))
  } else  {
    mcols_slots <- lapply(objects, function(x) x@elementMetadata)
    ## Will fail if not all the GAlignmentPairs objects in 'objects' have
    ## exactly the same metadata cols.
    ans_mcols <- do.call(rbind, mcols_slots)
  }

  ## Make 'ans' and return it.
  new(Class, NAMES=ans_NAMES,
      first=ans_first,
      last=ans_last,
      elementMetadata=ans_mcols)
}

setMethod("c", "GRangePairs",
          function(x, ..., ignore.mcols=FALSE, recursive=FALSE)
          {
            if (!identical(recursive, FALSE))
              stop("\"c\" method for GRangePairs objects ",
                   "does not support the 'recursive' argument")
            if (missing(x)) {
              objects <- list(...)
              x <- objects[[1L]]
            } else {
              objects <- list(x, ...)
            }
            combine_GRangePairs_objects(class(x), objects,
                                        use.names=FALSE,
                                        ignore.mcols=ignore.mcols)
          }
)

### -----------------------------------------------------------------
### swap method for GRangePairs: first becomes last and last becomes first
### Exported!
setGeneric("swap", function(x) standardGeneric("swap"))
setMethod("swap", "GRangePairs", function(x){
  BiocGenerics:::replaceSlots(x, first=last(x), last=first(x))
})

### -----------------------------------------------------------------
### unique: keep the unique GRangePairs
### Exported!
setMethod("unique", "GRangePairs", function(x){
  duplicatedIndex <- duplicated(paste(paste(first(x)), paste(last(x))))
  x[!duplicatedIndex]
})

### -----------------------------------------------------------------
### show methods for GRangePairs
###
.makeNakedMatFromGRangePairs <- function(x)
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  pair_cols_First <- cbind(seqnames=as.character(seqnames(x@first)),
                           strand=as.character(strand(x@first))
  )
  pair_cols_Last <- cbind(seqnames=as.character(seqnames(x@last)),
                          strand=as.character(strand(x@last))
  )
  x_first <- x@first
  first_cols <- cbind(ranges=showAsCell(ranges(x_first)))
  x_last <- x@last
  last_cols <- cbind(ranges=showAsCell(ranges(x_last)))
  ans <- cbind(pair_cols_First,
               #`:`=rep.int(":", lx),
               first_cols,
               `--`=rep.int("--", lx),
               pair_cols_Last,
               last_cols)
  if (nc > 0L) {
    tmp <- do.call(data.frame, lapply(mcols(x), showAsCell))
    ans <- cbind(ans, `|`=rep.int("|", lx), as.matrix(tmp))
  }
  ans
}

showGRangePairs <- function(x, margin="",
                            print.classinfo=FALSE,
                            print.seqinfo=FALSE){
  lx <- length(x)
  nc <- ncol(mcols(x))
  cat(class(x), " object with ",
      lx, " ", ifelse(lx == 1L, "pair", "pairs"),
      ", and ",
      nc, " metadata ", ifelse(nc == 1L, "column", "columns"),
      ":\n", sep="")
  out <- S4Vectors:::makePrettyMatrixForCompactPrinting(x,
                                                        .makeNakedMatFromGRangePairs)
  if(print.classinfo){
    .PAIR_COL2CLASS <- c(seqnames="Rle", strand="Rle")
    .HALVES_COL2CLASS <- c(ranges="IRanges")
    .COL2CLASS <- c(.PAIR_COL2CLASS, .HALVES_COL2CLASS, "--", 
                    .PAIR_COL2CLASS, .HALVES_COL2CLASS)
    classinfo <-
      S4Vectors:::makeClassinfoRowForCompactPrinting(x, .COL2CLASS)
    out <- rbind(classinfo, out)
  }
  if(nrow(out) != 0L)
    rownames(out) <- paste0(margin, rownames(out))
  ## We set 'max' to 'length(out)' to avoid the getOption("max.print")
  ## limit that would typically be reached when 'showHeadLines' global
  ## option is set to Inf.
  print(out, quote=FALSE, right=TRUE, max=length(out))
  if (print.seqinfo) {
    cat(margin, "-------\n", sep="")
    cat(margin, "seqinfo First: ", summary(seqinfo(x@first)), "\n", sep="")
    cat(margin, "seqinfo Last: ", summary(seqinfo(x@last)), "\n", sep="")
  }
}

setMethod("show", "GRangePairs",
          function(object)
            showGRangePairs(object, margin="  ",
                            print.classinfo=TRUE, print.seqinfo=TRUE)
)
