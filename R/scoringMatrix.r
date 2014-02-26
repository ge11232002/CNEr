### =================================================================
### pairwise genome alignment stuff
### -----------------------------------------------------------------

### -----------------------------------------------------------------
### Scoring matrix used in lastz from http://genomewiki.ucsc.edu/index.php/Hg19_100way_conservation_lastz_parameters.
### Exported!

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
  return(scoringMat)
}

### -----------------------------------------------------------------
### calculate the threshold score based on window size, identity and 
### scoring matrix
### Exported!
## This function tries to be compatible with the minScore_winSize thresholds
## used in previous CNEs identification.
## This function will output the equivilent threshold with a scoring matrix.
EScore <- function(scoringMat, winSize, minScore, gapOpen=NULL, gapExt=NULL){


}

