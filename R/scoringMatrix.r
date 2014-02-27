### =================================================================
### pairwise genome alignment stuff
### -----------------------------------------------------------------

### -----------------------------------------------------------------
### Scoring matrix used in lastz from http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_lastz_parameters.
### NOT Exported yet!

## This is the scoring matrix for highly divergent species.
HOXD55_MATRIX <- matrix(c(91, -90, -25, -100,
                          -90, 100, -100, -25,
                          -25, -100, 100, -90,
                          -100, -25, -90, 91),
                        nrow=4, ncol=4,
                        dimnames=list(c("A", "C", "G", "T"),
                                      c("A", "C", "G", "T"))
                        )

## This is the default medium scoring matrix.
HOXD70_MATRIX <- matrix(c(91, -114, -31, -123,
                          -114, 100, -125, -31,
                          -31, -125, 100,-114,
                          -123, -31, -114, 91),
                        nrow=4, ncol=4,
                        dimnames=list(c("A", "C", "G", "T"),
                                      c("A", "C", "G", "T"))
                        )

## This is the scoring matrix for close species.
HVSC_MATRIX <- matrix(c(90, -330, -236, -356,
                        -330, 100, -318, -236,
                        -236, -318, 100, -330,
                        -356, -236, -330, 90),
                      nrow=4, ncol=4,
                      dimnames=list(c("A", "C", "G", "T"),
                                    c("A", "C", "G", "T"))
                      )

### -----------------------------------------------------------------
### scoring matrix validator
### Not Exported!
normarg_scoringMat <- function(scoringMat){
  if(!is.matrix(scoringMat) || !is.numeric(scoringMat))
    stop("'", deparse(substitute(scoringMat)), "' must be a numeric matrix")
  if(!identical(rownames(scoringMat), DNA_BASES) || 
     !identical(colnames(scoringMat), DNA_BASES))
    stop("' rownames(", deparse(substitute(scoringMat)), ")' and 'colnames(", 
         deparse(substitute(scoringMat)),
         ")' must be the 4 DNA bases ('DNA_BASES')")
  if(any(is.na(scoringMat)))
    stop("'", deparse(substitute(scoringMat)), "' contains NAs")
  if(!isSymmetric(scoringMat))
    stop("'", deparse(substitute(scoringMat)), "' must be symmetric")
  return(scoringMat)
}

### -----------------------------------------------------------------
### calculate the threshold score based on window size, identity and 
### scoring matrix
### NOT Exported yet!
## This function tries to be compatible with the minScore_winSize thresholds
## used in previous CNEs identification.
## This function will output the equivilent threshold with a scoring matrix.
EScore <- function(scoringMat, winSize, minScore, 
                   ## NOT FINISHED!
                   bg=c(A=0.25, C=0.25, G=0.25, T=0.25),
                   gapOpen=NULL, gapExt=NULL){
  bg <- normargPriorParams(bg)
  scoringMat <- normarg_scoringMat(scoringMat)
  ans <- 0
  ## first calculate the scores from matches
  for(base in names(bg)){
    ans <- ans + minScore * bg[base] * scoringMat[base, base]
  }
  ## second calculate the scores from 

}

### Typical 'prior.params' vector: c(A=0.25, C=0.25, G=0.25, T=0.25)
### This is taken from Biostrings package.
normargPriorParams <- function(prior.params)
{
    if (!is.numeric(prior.params))
        stop("'prior.params' must be a numeric vector")
    if (length(prior.params) != length(DNA_BASES) ||
        !setequal(names(prior.params), DNA_BASES))
        stop("'prior.params' elements must be named A, C, G and T")
    ## Re-order the elements.
    prior.params <- prior.params[DNA_BASES]
    if (any(is.na(prior.params)) || any(prior.params < 0))
        stop("'prior.params' contains NAs and/or negative values")
    prior.params
}

