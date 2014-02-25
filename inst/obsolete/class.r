library(Biostrings)
library(GenomicRanges)

setClass(Class="axt",
         representation(targetNames="Rle",
                        targetRanges="IRanges",
                        targetSeqs="DNAStringSet",
                        queryNames="Rle",
                        queryRanges="IRanges",
                        querySeqs="DNAStringSet",
                        strand="Rle",
                        score="Rle"
                        ),
         prototype(
                   targetNames=Rle(factor()),
                   queryNames=Rle(factor()),
                   strand=Rle(strand())
                   )
         )

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Slot getters and setters.
###
setGeneric("targetNames", function(x) standardGeneric("targetNames"))
setMethod("targetNames", "axt", function(x) x@targetNames)
setGeneric("targetRanges", function(x) standardGeneric("targetRanges"))
setMethod("targetRanges", "axt", function(x) x@targetRanges)
setGeneric("targetSeqs", function(x) standardGeneric("targetSeqs"))
setMethod("targetSeqs", "axt", function(x) x@targetSeqs)
setGeneric("queryNames", function(x) standardGeneric("queryNames"))
setMethod("queryNames", "axt", function(x) x@queryNames)
setGeneric("queryRanges", function(x) standardGeneric("queryRanges"))
setMethod("queryRanges", "axt", function(x) x@queryRanges)
setGeneric("querySeqs", function(x) standardGeneric("querySeqs"))
setMethod("querySeqs", "axt", function(x) x@querySeqs)
setMethod("strand", "axt", function(x) x@strand)
setMethod("score", "axt", function(x) x@score)
setMethod("length", "axt", function(x) length(targetNames(x)))
setMethod("names", "axt", function(x) names(targetRanges(x)))

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Constructor.
###
axt = function(targetNames=Rle(), targetRanges=IRanges(),
               targetSeqs=DNAStringSet(),
               queryNames=Rle(), queryRanges=IRanges(),
               querySeqs=DNAStringSet(),
               strand=Rle("*", length(seqnames)), 
               score=Rle()){
  if(!is(targetNames, "Rle"))
    targetNames = Rle(targetNames)
  if(!is.factor(runValue(targetNames)))
    runValue(targetNames) = factor(runValue(targetNames))
  if(!is(queryNames, "Rle"))
    queryNames = Rle(queryNames)
  if(!is.factor(runValue(queryNames)))
    runValue(queryNames) = factor(runValue(queryNames))
  if(class(targetRanges) != "IRanges")
    targetRanges = as(targetRanges, "IRanges")
  if(class(queryRanges) != "IRanges")
    queryRanges = as(queryRanges, "IRanges")
  if(class(targetSeqs) != "DNAStringSet")
    targetSeqs = DNAStringSet(targetSeqs)
  if(class(querySeqs) != "DNAStringSet")
    querySeqs = DNAStringSet(querySeqs)
  if(!is(strand, "Rle"))
    strand = Rle(strand)
  if(!is.factor(runValue(strand)) || !identical(levels(runValue(strand)), levels(strand())))
    runValue(strand) = strand(runValue(strand))
  if(IRanges:::anyMissing(runValue(strand))){
    warning("missing values in strand converted to \"*\"")
    runValue(strand)[is.na(runValue(strand))] = "*"
  }
  if(!is(score, "Rle"))
    score = Rle(score)
  lx = max(length(targetNames), length(queryNames), 
            length(targetRanges), length(queryRanges),
            length(targetSeqs), length(querySeqs),
            length(strand), length(score))
  if(lx > 1){
    if(length(targetNames) == 1)
      targetNames = rep(targetNames, lx)
    if(length(queryNames) == 1)
      queryNames = rep(queryNames, lx)
    if(length(targetRanges) == 1)
      targetRanges = rep(targetRanges, lx)
    if(length(queryRanges) == 1)
      queryRanges = rep(queryRanges, lx)
    if(length(targetSeqs) == 1)
      targetSeqs = rep(targetSeqs, lx)
    if(length(querySeqs) == 1)
      querySeqs = rep(querySeqs, lx)
    if(length(strand) == 1)
      strand = rep(strand, lx)
    if(length(score) == 1)
      score = rep(score, lx)
  }
  new("axt", targetNames=targetNames, targetRanges=targetRanges, 
      targetSeqs=targetSeqs,
      queryNames=queryNames, queryRanges=queryRanges,
      querySeqs=querySeqs, strand=strand, score=score)
}

axt2 = function(target=GRanges(), query=GRanges(), score=integer()){
  new("axt2", target=target, query=query, score=score)
}

### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### Updating and cloning.
###
### An object is either 'update'd in place (usually with a replacement
### method) or 'clone'd (copied), with specified slots/fields overridden.

### For an object with a pure S4 slot representation, these both map to
### initialize. Reference classes will want to override 'update'. Other
### external representations need further customization.

setMethod("update", "axt",  # not exported
    function(object, ..., check=TRUE)
    {
        initialize(object, ...)
    }
)

setGeneric("clone", function(x, ...) standardGeneric("clone"))  # not exported
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
setMethod("[", "axt",
    function(x, i, ..., drop)
    {
        if(length(list(...)) > 0L)
          stop("invalid subsetting")
        if(missing(i))
          return(x)
        i <- IRanges:::normalizeSingleBracketSubscript(i, x)
        ans_targetNames = targetNames(x)[i]
        ans_targetRanges = targetRanges(x)[i]
        ans_targetSeqs = targetSeqs(x)[i]
        ans_queryNames = queryNames(x)[i]
        ans_queryRanges = queryRanges(x)[i]
        ans_querySeqs = querySeqs(x)[i]
        ans_strand = strand(x)[i]
        ans_score = score(x)[i]
        clone(x, targetNames=ans_targetNames,
                 targetRanges=ans_targetRanges,
                 targetSeqs=ans_targetSeqs,
                 queryNames=ans_queryNames,
                 queryRanges=ans_queryRanges,
                 querySeqs=ans_querySeqs,
                 strand=ans_strand,
                 score=ans_score)
    }
)

setMethod("[", "axt2",
          function(x, i, ..., drop){
            if(length(list(...)) > 0L)
              stop("invalid subsetting")
            if(missing(i))
              return(x)
            i = IRanges:::normalizeSingleBracketSubscript(i, x)
            ans_target = target(x)[i]
            ans_query = query(x)[i]
            ans_score = score(x)[i]
            clone(x, target=ans_target, query=ans_query, score=ans_score)
          }
          )
### - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -
### "show" method.
###
.makeNakedMatFromAxt = function(x){
  lx = length(x)
  ans = cbind(targetNames=as.character(targetNames(x)),
              targetRanges=IRanges:::showAsCell(targetRanges(x)),
              queryNames=as.character(queryNames(x)),
              queryRanges=IRanges:::showAsCell(queryRanges(x)),
              strand=as.character(strand(x)),
              score=as.character(score(x)),
              NULL)
  return(ans)
}

.rownames = function(names=NULL, len=NULL, tindex=NULL, bindex=NULL)
{
    if (is.null(tindex) & is.null(bindex)) {
        ## all lines
        if (len == 0L)
            character(0)
        else if (is.null(names))
            paste0("[", seq_len(len), "]")
        else
            names
    } else {
        ## head and tail
        if (!is.null(names)) {
            c(names[tindex], "...", names[bindex])
        } else {
            s1 = paste0("[", tindex, "]")
            s2 = paste0("[", bindex, "]")
            if (all(tindex == 0))
                s1 = character(0)
            if (all(bindex == 0))
                s2 = character(0)
            c(s1, "...", s2)
        }
    }
}

makePrettyMatrixForCompactPrintingAxt = function(x, makeNakedMat.FUN){
  nhalf = 5L
  lx = length(x)
  if(is.null(nhead <- getOption("showHeadLines")))
    nhead = nhalf
  if(is.null(ntail <- getOption("showTailLines")))
    ntail = nhalf
  if(lx < (nhead+ntail+1L)){
    ans = makeNakedMat.FUN(x)
    ans_rownames = .rownames(names(x), lx)
  }else{
    top_idx = 1:nhead
    if(nhead == 0)
      top_idx = 0
    bottom_idx=(lx-ntail+1L):lx
    if(ntail == 0)
      bottom_idx = 0
    ans_top = makeNakedMat.FUN(x[top_idx])
    ans_bottom = makeNakedMat.FUN(x[bottom_idx])
    ans = rbind(ans_top,
                 matrix(rep.int("...", ncol(ans_top)), nrow=1L),
                 ans_bottom)
    ans_rownames = .rownames(names(x), lx, top_idx, bottom_idx)
  }
  rownames(ans) = format(ans_rownames, justify="right")
  return(ans)
}

showAxt = function(x, margin=""){
  lx = length(x)
  cat(class(x), " with ",
      lx, " ", ifelse(lx == 1L, "alignment pair", "alignment pairs"),
      ":\n", sep="")
  out = makePrettyMatrixForCompactPrintingAxt(x, .makeNakedMatFromAxt)
  if(nrow(out) != 0L)
    rownames(out) = paste0(margin, rownames(out))
  print(out, quote=FALSE, right=TRUE)
}

setMethod("show", "axt",
          function(object)
            showAxt(object, margin="  ")
)






