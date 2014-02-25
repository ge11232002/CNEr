#include "R.h"
#include  <ctype.h>  //for the tolower()
#include "ucsc/common.h"
#include "ucsc/linefile.h"
#include "ucsc/hash.h"
#include "ucsc/obscure.h"
#include "ucsc/options.h"
#include "ucsc/axt.h"
#include "Rdefines.h"
#include "Rinternals.h"
#include <string.h>
#include "IRanges_interface.h"
#include "ucsc/binRange.h"
#include "XVector_interface.h"
#include <R_ext/Rdynload.h>

/* Data structures to represent start and end coordinate pairs.
 *  * Used to store filters in memory. */
struct range
/* Start and end coordinate pair */
{
  int start;    /* Start (0 based) */
  int end;    /* End (non-inclusive) */
};

struct rangeArray
/* Array of start and end coordinate pairs */
{
  int n;
  struct range *ranges;
};

struct slRange
/* Start and end coordinate pair as linked list item */
{
  struct slRange *next;
  int start;    /* Start (0 based) */
  int end;    /* End (non-inclusive) */
};


/* Data structure used to represent different thresholds and 
 * intermediate results for each */

struct slThreshold
{
  struct slThreshold *next;
  int minScore;
  int winSize;
  int ceStart;
  int ceEnd;
  int nrCNE;
  struct slCNE *CNE;
  FILE *outFile;
};

struct slCNE
{
  struct slCNE *next;
  char *tName; // Name of the target sequence.
  int tStart; // The 1-based coordinate
  int tEnd; // The 1-based coordinate
  char *qName; // Name of the query sequence.
  int qStart;
  int qEnd;
  char strand;
  float score;
  char *cigar;
};

struct slAllCNE
{
  struct slAllCNE *next;
  int minScore;
  int winSize;
  struct slCNE *CNE;
};

/* io.c */
SEXP myReadBed(SEXP filepath);

SEXP myReadAxt(SEXP filepath);

SEXP axt_info(SEXP filepath);

/* utils.c */
SEXP bin_from_coord_range(SEXP starts, SEXP ends);

SEXP bin_ranges_from_coord_range(SEXP start, SEXP end);

/* ceScan.c */
SEXP myCeScan(SEXP tFilterNames, SEXP tFilterStarts, SEXP tFilterEnds,
              SEXP qFilterNames, SEXP qFilterStarts, SEXP qFilterEnds,
              SEXP sizeNames, SEXP sizeSizes, SEXP axttNames,
              SEXP axttStart, SEXP axttEnd, SEXP axttStrand,
              SEXP axttSym, SEXP axtqNames, SEXP axtqStart,
              SEXP axtqEnd, SEXP axtqStrand, SEXP axtqSym,
              SEXP score, SEXP symCount, SEXP winSize,
              SEXP minScore, SEXP outputFiles);

SEXP myCeScanFile(SEXP axtFiles, SEXP tFilterFile, SEXP qFilterFile,
                SEXP sizeNames, SEXP sizeSizes,
                SEXP winSize, SEXP minScore, SEXP outputFiles);

