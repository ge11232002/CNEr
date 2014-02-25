
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
                 )
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
              if(any(object@symCount <= 0L))
                return("Then symCount must be larger than 0!")
              return(TRUE)
            }
            )

### -----------------------------------------------------------------
### CNE class
### 
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
###
setMethod("targetRanges", "Axt", function(x) x@targetRanges)
setMethod("targetSeqs", "Axt", function(x) x@targetSeqs)
setMethod("queryRanges", "Axt", function(x) x@queryRanges)
setMethod("querySeqs", "Axt", function(x) x@querySeqs)
setMethod("score", "Axt", function(x) x@score)
setMethod("symCount", "Axt", function(x) x@symCount)
setMethod("length", "Axt", function(x) length(targetRanges(x)))

### -----------------------------------------------------------------
### CNE Slot getters and setters.
###
setMethod("assembly1", "CNE", function(x) x@assembly1)
setMethod("assembly2", "CNE", function(x) x@assembly2)
setMethod("CNE1", "CNE", function(x) x@CNE1)
setMethod("CNE2", "CNE", function(x) x@CNE2)
setMethod("thresholds", "CNE", function(x) x@thresholds)
setMethod("CNEMerged", "CNE", function(x) x@CNEMerged)
setMethod("CNERepeatsFiltered", "CNE", function(x) x@CNERepeatsFiltered)

### -- - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Axt Constructor.
###
Axt <- function(targetRanges=GRanges(), targetSeqs=DNAStringSet(),
               queryRanges=GRanges(), querySeqs=DNAStringSet(),
               score=integer(), symCount=integer()){
  new("Axt", targetRanges=targetRanges, targetSeqs=targetSeqs,
      queryRanges=queryRanges, querySeqs=querySeqs,
      score=score, symCount=symCount)
}

### -----------------------------------------------------------------
### CNE constructor.
###
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

setMethod("update", "Axt",
          function(object, ..., check=TRUE){
            initialize(object, ...)
          }
          )

setMethod("update", "CNE",
          function(object, ..., check=TRUE){
            initialize(object, ...)
          }
          )

setMethod("clone", "ANY",  # not exported
    function(x, ...)
    {
        if (nargs() > 1L)
            initialize(x, ...)
        else
            x
    }
)


### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Subsetting and combining.
###

setMethod("[", "Axt",
          function(x, i, ..., drop){
            if(length(list(...)) > 0L)
              stop("invalid subsetting")
            if(missing(i))
              return(x)
            #i = IRanges:::normalizeSingleBracketSubscript(i, x)
            ans_targetRanges <- targetRanges(x)[i]
            ans_targetSeqs <- targetSeqs(x)[i]
            ans_queryRanges <- queryRanges(x)[i]
            ans_querySeqs <- querySeqs(x)[i]
            ans_score <- score(x)[i]
            ans_symCount <- symCount(x)[i]
            clone(x, targetRanges=ans_targetRanges, targetSeqs=ans_targetSeqs,
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


setMethod("subAxt", "Axt",
## This is to fetch the Axts within the specific chrs, starts, ends 
## based on target sequences.
          function(x, chr, start, end, #strand=c("+", "-", "*"), 
                   select=c("target", "query"),
                   type=c("any", "within"),
                   qSize=NULL){
            type <- match.arg(type)
            select <- match.arg(select)
            if(select == "query"){
              if(is.null(qSize))
                stop("When selecting on query alignments, 
                     qSize must be provided.")
              if(!is(qSize, "integer"))
                stop("qSize must be an integer object.")
            }
            if(select == "target"){
              searchGRanges = GRanges(seqnames=chr,
                                      ranges=IRanges(start=start, end=end),
                                      strand="+")
              if(length(searchGRanges) == 0)
                return(x)
              index <- queryHits(findOverlaps(targetRanges(x), 
                                             searchGRanges, type=type, 
                                             select="all"))
            }else if(select == "query"){
              # first search Axts on positive strand
              searchGRanges = GRanges(seqnames=chr,
                                      ranges=IRanges(start=start, end=end),
                                      strand="+")
              if(length(searchGRanges) == 0)
                return(x)
              indexPositive <- queryHits(findOverlaps(queryRanges(x),
                                                     searchGRanges, type=type,
                                                     select="all"))
              # then search Axts on negative strand.
              # we need to prepare the searchGRanges on negative strand.
              searchGRanges <- GRanges(seqnames=chr,
                                      ranges=IRanges(start=qSize-end+1,
                                                     end=qSize-start+1),
                                      strand="-")
              indexNegative <- queryHits(findOverlaps(queryRanges(x),
                                                     searchGRanges, type=type,
                                                     select="all"))
              index <- sort(c(indexPositive, indexNegative))
            }else{
              stop("Wrong select!")
            }
            return(x[index])
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

