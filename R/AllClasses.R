
### -----------------------------------------------------------------
### Axt class
### Exported!
setClass(Class="Axt",
         slots=c(targetRanges="GRanges",
                 targetSeqs="DNAStringSet",
                 queryRanges="GRanges",
                 querySeqs="DNAStringSet",
                 score="integer",
                 symCount="integer"
                 ),
         prototype=list(targetRanges=GRanges(),
                        targetSeqs=DNAStringSet(),
                        queryRanges=GRanges(),
                        querySeqs=DNAStringSet(),
                        score=integer(0),
                        symCount=integer(0))
         )

setValidity("Axt",
            function(object){
              if(!isConstant(c(length(object@targetRanges), 
                             length(object@targetSeqs),
                             length(object@queryRanges), 
                             length(object@querySeqs),
                             length(object@score), 
                             length(object@symCount))))
                return("The lengths of targetRanges, targetSeqs,
                       queryRanges, querySeqs, score and symCount
                       must have be same!")
              if(any(object@symCount < 0L))
                return("Then symCount must be equal or larger than 0!")
              ## Test the class
              if(class(object@targetRanges) != "GRanges")
                return("'x@targetRanges' must be a GRanges instance")
              if(class(object@queryRanges) != "GRanges")
                return("'x@queryRanges' must be a GRanges instance")
              if(class(object@targetSeqs) != "DNAStringSet")
                return("'x@targetSeqs' must be a DNAStringSet instance")
              if(class(object@querySeqs) != "DNAStringSet")
                return("'x@querySeqs' must be a DNAStringSet instance")
              return(TRUE)
            }
            )

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

### -----------------------------------------------------------------
### CNE class
### Exported!
setClass(Class="CNE",
         slots=c(assembly1="character",
                 assembly2="character",
                 thresholds="character",
                 CNE1="list",
                 CNE2="list",
                 CNEMerged="list",
                 CNERepeatsFiltered="list",
                 alignMethod="character"
                 )
         )

setValidity("CNE",
            function(object){
              if(length(object@assembly1) != 1L)
                return("The name of assembly1 must be length 1!")
              if(length(object@alignMethod) != 1L)
                return("The align method must be length 1!")
              if(length(object@assembly2) != 1L)
                return("The name of assembly2 must be length 1!")
              if(!all(grepl("^\\d+_\\d+$", object@thresholds)))
                return("The thresholds must be in format of 49_50!")
              if(any(as.integer(
                       sapply(strsplit(object@thresholds, "_"), "[", 2))
                 < as.integer(
                       sapply(strsplit(object@thresholds, "_"), "[", 1))))
                return("The window size cannot be smaller than identity score!")
              if(length(object@CNE1) != length(object@thresholds) ||
                 length(object@CNE2) != length(object@thresholds) ||
                 length(object@CNEMerged) != length(object@thresholds) ||
                 length(object@CNERepeatsFiltered) != length(object@thresholds))
                return("The number of cne tables must be same with
                       number of thresholds!")
              return(TRUE)
            }
            )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Axt Slot getters and setters.
### Exported!
setMethod("targetRanges", "Axt", function(x) x@targetRanges)
setMethod("targetSeqs", "Axt", function(x) x@targetSeqs)
setMethod("queryRanges", "Axt", function(x) x@queryRanges)
setMethod("querySeqs", "Axt", function(x) x@querySeqs)
setMethod("score", "Axt", function(x) x@score)
setMethod("symCount", "Axt", function(x) x@symCount)
setMethod("nchar", "Axt", function(x) x@symCount)
setMethod("length", "Axt", function(x) length(targetRanges(x)))

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
### CNE Slot getters and setters.
### Exported!
setMethod("assembly1", "CNE", function(x) x@assembly1)
setMethod("assembly2", "CNE", function(x) x@assembly2)
setMethod("CNE1", "CNE", function(x) x@CNE1)
setMethod("CNE2", "CNE", function(x) x@CNE2)
setMethod("thresholds", "CNE", function(x) x@thresholds)
setMethod("CNEMerged", "CNE", function(x) x@CNEMerged)
setMethod("CNERepeatsFiltered", "CNE", function(x) x@CNERepeatsFiltered)

### -----------------------------------------------------------------
### Axt Constructor.
### Exported!
Axt <- function(targetRanges=GRanges(), targetSeqs=DNAStringSet(),
               queryRanges=GRanges(), querySeqs=DNAStringSet(),
               score=integer(0), symCount=integer(0)){
  new("Axt", targetRanges=targetRanges, targetSeqs=targetSeqs,
      queryRanges=queryRanges, querySeqs=querySeqs,
      score=score, symCount=symCount)
}

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
### CNE constructor.
### Exported!
CNE <- function(assembly1=character(), assembly2=character(),
                thresholds=character(),
                CNE1=list(), CNE2=list(),
                CNEMerged=list(), CNERepeatsFiltered=list(),
                alignMethod=character()
                ){
  new("CNE", assembly1=assembly1, assembly2=assembly2,
      thresholds=thresholds, CNE1=CNE1, CNE2=CNE2,
      CNEMerged=CNEMerged, CNERepeatsFiltered=CNERepeatsFiltered,
      alignMethod=alignMethod)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.

### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

# setMethod("update", "Axt",
#           function(object, ..., check=TRUE){
#             initialize(object, ...)
#           }
#           )
# 
# setMethod("update", "CNE",
#           function(object, ..., check=TRUE){
#             initialize(object, ...)
#           }
#           )
# 
# setMethod("clone", "ANY",  # not exported
#     function(x, ...)
#     {
#         if (nargs() > 1L)
#             initialize(x, ...)
#         else
#             x
#     }
# )


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
###

setMethod("[", "Axt",
          function(x, i, ..., drop){
            if(length(list(...)) > 0L)
              stop("invalid subsetting")
            if(missing(i))
              return(x)
            i <- normalizeSingleBracketSubscript(i, x)
            ans_targetRanges <- targetRanges(x)[i]
            ans_targetSeqs <- targetSeqs(x)[i]
            ans_queryRanges <- queryRanges(x)[i]
            ans_querySeqs <- querySeqs(x)[i]
            ans_score <- score(x)[i]
            ans_symCount <- symCount(x)[i]
            initialize(x, targetRanges=ans_targetRanges, targetSeqs=ans_targetSeqs,
                  queryRanges=ans_queryRanges, querySeqs=ans_querySeqs,
                  score=ans_score, symCount=ans_symCount)
          }
          )

setMethod("c", "Axt",
          function(x, ...){
            if(missing(x)){
              args <- unname(list(...))
              x <- args[[1L]]
            }else{
              args <- unname(list(x, ...))
            }
            if(length(args) == 1L)
              return(x)
            arg_is_null <- sapply(args, is.null)
            if(any(arg_is_null))
              args[arg_is_null] <- NULL
            if(!all(sapply(args, is, class(x))))
              stop("all arguments in '...' must be ", 
                   class(x), " objects (or NULLs)")
            new_targetRanges <- do.call(c, lapply(args, targetRanges))
            new_targetSeqs <- do.call(c, lapply(args, targetSeqs))
            new_queryRanges <- do.call(c, lapply(args, queryRanges))
            new_querySeqs <- do.call(c, lapply(args, querySeqs))
            new_score <- do.call(c, lapply(args, score))
            new_symCount <- do.call(c, lapply(args, symCount))

            initialize(x,
                       targetRanges=new_targetRanges,
                       targetSeqs=new_targetSeqs,
                       queryRanges=new_queryRanges,
                       querySeqs=new_querySeqs,
                       score=new_score,
                       symCount=new_symCount)
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

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###
### 'x' must be an XString or MaskedXString object.
toSeqSnippet <- function(x, width)
{
    if (width < 7L)
        width <- 7L
    seqlen <- length(x)
    if (seqlen <= width) {
        as.character(x)
    } else {
        w1 <- (width - 2) %/% 2
        w2 <- (width - 3) %/% 2
        paste(as.character(subseq(x, start=1, width=w1)),
              "...",
              as.character(subseq(x, end=seqlen, width=w2)),
              sep="")
    }
}


.axt.show_frame_line <- function(x, i, iW, tNameW, tStartW, tEndW, 
                                qNameW, qStartW, qEndW, scoreW){
  cat(format(i, width=iW, justify="right"), " ",
      format(as.character(seqnames(targetRanges(x)[i])), 
             width=tNameW, justify="right"), " ",
      format(start(targetRanges(x)[i]), width=tStartW, justify="right"), " ",
      format(end(targetRanges(x)[i]), width=tEndW, justify="right"), " ",
      format(as.character(seqnames(queryRanges(x)[i])), 
             width=qNameW, justify="right"), " ",
      format(start(queryRanges(x)[i]), width=qStartW, justify="right"), " ",
      format(end(queryRanges(x)[i]), width=qEndW, justify="right"), " ",
      format(as.character(strand(queryRanges(x))[i]), 
             width=1, justify="right"), " ",
      format(score(x)[i], width=scoreW, justify="right"), " ",
      sep=""
      )
  cat("\n")
  snippetWidth <- getOption("width")
  seq_snippet <- toSeqSnippet(targetSeqs(x)[[i]], snippetWidth)
  cat(seq_snippet)
  cat("\n")
  seq_snippet <- toSeqSnippet(querySeqs(x)[[i]], snippetWidth)
  cat(seq_snippet)
  cat("\n")
}

showAxt <- function(x, margin="", half_nrow=5L){
  lx <- length(x)
  if(is.null((head_nrow = getOption("showHeadLines"))))
    head_nrow = half_nrow
  if(is.null((tail_nrow = getOption("showTailLines"))))
    tail_nrow = half_nrow
  iW = nchar(as.character(lx))
  if(lx < (2*half_nrow+1L) | (lx < (head_nrow+tail_nrow+1L))) {
    tNameW <- max(nchar(as.character(seqnames(targetRanges(x)))))
    tStartW <- max(nchar(as.character(start(targetRanges(x)))))
    tEndW <- max(nchar(as.character(end(targetRanges(x)))))
    qNameW <- max(nchar(as.character(seqnames(queryRanges(x)))))
    qStartW <- max(nchar(as.character(start(queryRanges(x)))))
    qEndW <- max(nchar(as.character(end(queryRanges(x)))))
    scoreW <- max(nchar(as.character(score(x))))
    for(i in seq_len(lx))
      .axt.show_frame_line(x, i, iW, tNameW, tStartW, tEndW, 
                           qNameW, qStartW, qEndW, scoreW)
  }else{
    tNameW <- max(nchar(as.character(seqnames(targetRanges(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    tStartW <- max(nchar(as.character(start(targetRanges(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    tEndW <- max(nchar(as.character(end(targetRanges(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    qNameW <- max(nchar(as.character(seqnames(queryRanges(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    qStartW <- max(nchar(as.character(start(queryRanges(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    qEndW <- max(nchar(as.character(end(queryRanges(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)]))))
    scoreW <- max(nchar(as.character(score(x)
                        [c(1:head_nrow, (lx-tail_nrow+1L):lx)])))
    if(head_nrow > 0){
      for(i in 1:head_nrow)
        .axt.show_frame_line(x, i, iW, tNameW, tStartW, tEndW, 
                             qNameW, qStartW, qEndW, scoreW)
    }
    cat(format("...", width=iW, justify="right"),
        format("...", width=tNameW, justify="right"),
        format("...", width=tStartW, justify="right"),
        format("...", width=tEndW, justify="right"),
        format("...", width=qNameW, justify="right"),
        format("...", width=qStartW, justify="right"),
        format("...", width=qEndW, justify="right"),
        format("...", width=scoreW, justify="right")
        )
    cat("\n")
    if(tail_nrow > 0){
      for(i in (lx-tail_nrow+1L):lx)
        .axt.show_frame_line(x, i, iW, tNameW, tStartW, tEndW, 
                             qNameW, qStartW, qEndW, scoreW)
    }
  }
}
    #out = makePrettyMatrixForCompactPrintingAxt(x, .makeNakedMatFromAxt)
    #if(nrow(out) != 0L)
    #      rownames(out) = paste0(margin, rownames(out))
    #  print(out, quote=FALSE, right=TRUE)

setMethod("show", "Axt",
          function(object){
            lx <- length(object)
            cat(" A ", class(object), " with ", length(object), " ", 
                ifelse(lx == 1L, "alignment pair", "alignment pairs"), 
                ":\n", sep="")
            if(lx != 0){
              showAxt(object, margin="  ")
            }
          }
)

### -----------------------------------------------------------------
### show methods for GRangePairs
###

.makeNakedMatFromGRangePairs <- function(x)
{
  lx <- length(x)
  nc <- ncol(mcols(x))
  pair_cols <- cbind(seqnames=as.character(seqnames(x)),
                     strand=as.character(strand(x)))
  x_first <- x@first
  first_cols <- cbind(ranges=extractROWS(ranges(x_first)))
  x_last <- x@last
  last_cols <- cbind(ranges=showAsCell(ranges(x_last)))
  ans <- cbind(pair_cols,
               `:`=rep.int(":", lx),
               first_cols,
               `--`=rep.int("--", lx),
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
  
}

setMethod("show", "GRangePairs",
          function(object)
            showGRangePairs(object, margin="  ",
                                print.classinfo=TRUE, print.seqinfo=TRUE)
)
