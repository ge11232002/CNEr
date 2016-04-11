### -----------------------------------------------------------------
### GRangePairs: this copies the implementation of GAlignmentPairs
### Exported!
setClass(Class="GRangePairs",
         contains="List",
         slots=c(NAMES="characterORNULL",      # R doesn't like @names !!
                 first="GRanges",
                 last="GRanges",
                 elementMetadata="DataFrame"),
         prototype=list(elementType="GRanges"))

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
          function(x, recursive=TRUE, use.names=TRUE)
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
