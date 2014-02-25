/* ceScan.c - scan axt alignment for conserved elements */

#include "common.h"
#include "linefile.h"
#include "hash.h"
#include "obscure.h"
#include "options.h"
#include "axt.h"


/********************************************
 *** DATA STRUCTURES AND GLOBAL VARIABLES ***
 ********************************************/

/* Commandline options */

static struct optionSpec options[] = {
   {"tFilter", OPTION_STRING},
   {"qFilter", OPTION_STRING},
   {"outPrefix", OPTION_STRING},
   {NULL, 0},
};

/* Scoring matrix.
 * This will be set by setBpScores() to 1 for matches and 0 for mismatches and gaps. */

#define NR_CHARS 128
typedef int bpScores_t[NR_CHARS][NR_CHARS];
static bpScores_t bpScores;

/* Data structures to represent start and end coordinate pairs. 
 * Used to store filters in memory. */

struct range 
/* Start and end coordinate pair */
{
  int start;		/* Start (0 based) */
  int end;		/* End (non-inclusive) */
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
  int start;		/* Start (0 based) */
  int end;		/* End (non-inclusive) */
};


/* Data structure used to represent different thresholds and intermediate results for each */

struct slThreshold
{
  struct slThreshold *next;
  int minScore;
  int winSize;
  int ceStart;
  int ceEnd;
  FILE *outFile;
};


/*****************
 *** FUNCTIONS ***
 *****************/

void usage()
/* Explain usage and exit. */
{
errAbort(
  "ceScan - Find conserved elements in axt alignments\n\n"
  "Usage:\n"
  "   ceScan [options] query-sizes-file thresholds axt-files\n\n"
  "The program detects conserved elements as regions with at least I identities\n"
  "over C alignment columns. Conserved elements that overlap in an alignment\n"
  "(i.e. overlap on both genomes) are merged. The ends of conserved elements are\n"
  "trimmed from mismatches and gaps. Reported conserved elements may therefore be\n"
  "shorter than C (if I < C). If filter files are given, the regions specified in\n"
  "those files will be excluded from the scan, so that reported conserved\n"
  "elements do not overlap filter regions.\n\n"
  "Required Arguments:\n"
  "   query-sizes-file - Sizes of the query sequences, as obtained from twoBitInfo.\n"
  "   thresholds - A threshold is specified as I,C (e.g. 45,50 for 45 (95%%)\n"
  "                identities over 50 columns).\n"
  "                Multiple thresholds can be given, separated by whitespace.\n"
  "                The program will look for conserved elements at each threshold\n"
  "                given. If you specify multiple thresholds, you should also\n"
  "                use the -outPrefix option, so that results for different\n"
  "                thresholds are written to different files.\n"
  "   axt-files - Alignment files to scan. Multiple files can be given.\n"
  "\n"
  "Options:\n"
  "   -tFilter=bedfile - Regions of target assembly to exclude from scan.\n"
  "   -qFilter=bedfile - Regions of query assembly to exclude from scan.\n"
  "                      The bed files can be unsorted and redundant.\n"
  "                      Only the first 3 columns will be read.\n"
  "   -outPrefix=prefix - Prefix for output filenames. This will be combined\n"
  "                       with threshold parameters to make filenames:\n"
  "                       prefix_minScore_winSize\n"
  "                       By default all output is written to stdout.\n"
  "\n"
  "Output format:\n"
  "   Tab-separated file(s) with one line for each conserved element.\n"
  "   Columns:\n"
  "     target sequence name, target start coordinate, target end coordinate,\n"
  "     query sequence name, query start coordinate, query end coordinate,\n"
  "     alignment orientation (+ or -), percent identity, CIGAR string.\n"
  "   The coordinate ranges are half-open, zero-based.\n"
  "   The percent identity is computed as percentage of alignment columns\n"
  "   with identities.\n"
  "   The CIGAR string format follows the DAS specification for alignments\n"
  "   (see http://www.dasregistry.org/extension_alignment.jsp).\n"
  "\n"
  );
}


struct hash *loadIntHash(char *fileName)
/* Read in a file full of name/number lines into a hash keyed
 * by name with number values. Adapted from axtToMaf.c. */
{
  struct lineFile *lf = lineFileOpen(fileName, TRUE);
  char *row[2];
  struct hash *hash = newHash(0);
  
  while (lineFileRow(lf, row)) {
    int num = lineFileNeedNum(lf, row, 1);
    hashAddInt(hash, row[0], num);
  }
  
  lineFileClose(&lf);
  return hash;
}


struct hash *readBed(char *fileName)
/* Read a 3-column bed file into a hash, where keys are sequence names
 * and values are linked lists of coordinate ranges (slRange structures). */
{
  struct lineFile *lf = lineFileOpen(fileName, TRUE);
  struct hash *hash = newHash(0);
  struct hashEl *hel;
  struct slRange *range;
  char *row[3];

  while (lineFileRow(lf, row)) {
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    AllocVar(range);
    range->next = NULL;
    range->start = lineFileNeedNum(lf, row, 1);
    range->end = lineFileNeedNum(lf, row, 2);
    if (range->start > range->end)
      errAbort("start after end line %d of %s", lf->lineIx, lf->fileName);
    
    hel = hashLookup(hash, row[0]);
    if (hel == NULL)
      hel = hashAdd(hash, row[0], range);
    else {
      slSafeAddHead(&hel->val, range);
    }
  }

  lineFileClose(&lf);

  return hash;
}


int slRangeCmpStart(const void *va, const void *vb)
/* Comparison function to sort linked list of ranges by start coordinate. */
{
  const struct slRange *a = *((struct slRange **)va);
  const struct slRange *b = *((struct slRange **)vb);
  return a->start - b->start;
}


void collapseRangeList(struct hashEl *hel)
/* Collapse a range list to a set of sorted, non-overlapping ranges. */
{
  struct slRange *a, *b;
  slSort(&hel->val, slRangeCmpStart); /* sort by start coord */
  a = hel->val;
  while((b = a->next)) {
    if(b->start <= a->end) {
      if(a->end < b->end) a->end = b->end;
      a->next = b->next;
      freez(&b);
    }
    else a = b;
  }
  /*for(a = hel->val; a; a=a->next) {
    printf("%d\t%d\n", a->start, a->end);
    }*/
}


void convertRangeListToArray(struct hashEl *hel)
/* Convert a linked list of ranges to an array. 
 * The reason for doing this is that we can do a fast binary search on the array. */
{
  struct slRange *list, *slEl;
  struct range *arrayEl;
  struct rangeArray *arrayInfo;
  int n;

  list = hel->val;
  n = slCount(list)+1;
  AllocVar(arrayInfo);
  arrayInfo->n = n;
  arrayInfo->ranges = arrayEl = needMem(n * sizeof(*arrayEl));
  hel->val = arrayInfo;

  while((slEl = slPopHead(&list))) {
    arrayEl->start = slEl->start;
    arrayEl->end = slEl->end;
    free(slEl);
    arrayEl++;
  }

  /* The last array element is a "dummy" element that contains a coordinate pair
   * beyond any chromosome size. The presence of this element simplifies going 
   * through the array in scanAxt() as it removes the need for an out-of-bounds check. */
  arrayEl->start = 1e9;
  arrayEl->end = 1e9+1;
}


void printRangeArray(struct hashEl *hel)
/* Print a range array. For debugging purposes only. */
{
  struct rangeArray *arrayInfo = hel->val;
  struct range *ranges = arrayInfo->ranges; 
  int i;
  printf("%s n=%d\n", hel->name, arrayInfo->n);
  for(i = 0; i < arrayInfo->n; i++) {
    printf("%02d: %d - %d\n", i, ranges[i].start, ranges[i].end);
  }
}


struct range *searchRangeArray(struct rangeArray *arrayInfo, int key)
/* Binary search range array. */
{
  struct range *array = arrayInfo->ranges;
  int low = 0; 
  int high = arrayInfo->n - 1;
  int mid;

  while(low <= high) {
    mid = (low+high)/2;
    if(key <= array[mid].start) high = mid - 1;
    else if(key > array[mid].end) low = mid + 1;
    else return array+mid; /* return pointer to range that contains key */
  }

  /* key not found: return pointer to nearest higher range or abort if there is no higher range 
   * (there should be one because we have added a dummy range with very high values) */
  if(low >= arrayInfo->n) errAbort("searchRangeArray: key %d out of bounds\n", key);
  return array+low;
}


struct hash *readFilter(char *fileName)
/* Load a filter file. */
{
  struct hash *hash = readBed(fileName);
  hashTraverseEls(hash, collapseRangeList);
  hashTraverseEls(hash, convertRangeListToArray);
  /* hashTraverseEls(hash, printRangeArray); */
  return hash;
}


struct hash *makeReversedFilter(struct hash *f1, struct hash *chrSizes)
/* Given a filter, create a reversed filter where coordinates increase in the opposite direction.
 * We use this for filtering alignments that have qStrand == '-'. */
{
  struct hash *f2 = newHash(0);
  struct hashCookie cookie = hashFirst(f1);
  struct hashEl *hel;
  struct rangeArray *fwd, *rev;
  struct range *arrayEl;
  int i, j, n, chrSize;

  /* Iterate over all sequences (chromosomes) in filter */
  while((hel = hashNext(&cookie))) {  

    /* get sequence size */
    chrSize = hashIntVal(chrSizes, hel->name);

    /* get forward range array */
    fwd = hel->val;

    /* allocate memory for reversed range array */
    AllocVar(rev);
    n = rev->n = fwd->n; /* set nr of elements in range */
    rev->ranges = arrayEl = needMem(n * sizeof(struct range));

    /* copy dummy range */
    rev->ranges[n-1] = fwd->ranges[n-1];

    /* reverse other ranges */
    for(i = 0, j = n-2; j >= 0; i++, j--) {
      rev->ranges[i].start = chrSize - fwd->ranges[j].end;
      rev->ranges[i].end = chrSize - fwd->ranges[j].start;
    }

    /* add range array to hash keyed by sequence name */
    hashAdd(f2, hel->name, rev);
  }

  /* return reverse filter */
  return f2;
} 


struct range *searchFilter(struct hash *filter, char *chrom, int pos)
/* Find the first filter at or following a given position */
{
  struct hashEl *hel;

  hel = hashLookup(filter, chrom);   /* find range array for sequence (chromosome) */
  if(hel) return searchRangeArray(hel->val, pos); /* search range array by position */
  else return NULL;
}


void setBpScores(bpScores_t ss)
/* Set scoring matrix to 1 for matches and 0 for mismatches and gaps. */
{
  unsigned int i, j;
  int a, A;
  char validChars[] = "ACGT";

  // printf("%d\n", (int) sizeof(bpScores_t));

  for (i = 0; i < NR_CHARS; ++i)
    for (j = 0; j < NR_CHARS; ++j)
      ss[i][j] = 0;
  for (i = 0; i < sizeof(validChars)-1; ++i) {
    A = validChars[i];
    a = tolower(A);
    ss[A][A] = ss[a][A] = ss[A][a] = ss[a][a] = 1;
  }
}


void printCigarString(FILE *fh, struct axt *axt, int i, int j)
/* Print CIGAR string that summarizes alignment */
{
  char type = 'M'; /* in our case first column is always match */
  char newType;
  int count = 0;

  for(; i <= j; i++) {
    /* Determine column type */
    if(axt->tSym[i] == '-') newType = 'D';
    else if(axt->qSym[i] == '-') newType = 'I';
    else newType = 'M';
    /* If same type as previous, just increase count, otherwise output previous */
    if(type == newType) count++;
    else {
      fprintf(fh, "%d%c", count, type);
      type = newType;
      count = 1;
    }   
  }
  
  if(count) fprintf(fh, "%d%c", count, type);
}


void printElement(struct slThreshold *tr, struct axt *axt, struct hash *qSizes, int *profile, int *tPosList, int *qPosList)
/* Print one conserved element on stdout.
 * Arguments:
 * tr - contains threshold-specific information:
 *      parameters used to find CE, CE location in alignment, and filehandle to print to
 * axt - alignment
 * qSizes - query assembly chromosome sizes
 * profile - cumulative conservation profile for alignment
 * tPosList, qPosList - target and query position arrays for alignment  
 */
{
  int score, qStart, qEnd, qSize;
  int i = tr->ceStart; /* start column of conserved element in alignment */
  int j = tr->ceEnd; /* end column of conserved element in alignment */
  
  /* truncate edges (mismatches and gaps) */
  while(bpScores[ (int) axt->qSym[i] ][ (int) axt->tSym[i] ] <= 0) i++;
  while(bpScores[ (int) axt->qSym[j] ][ (int) axt->tSym[j] ] <= 0) j--;

  /* compute score */
  score = profile[j] - profile[i] + bpScores[ (int) axt->qSym[i] ][ (int) axt->tSym[i] ];

  /* recompute query positions if query strand is - */
  if(axt->qStrand == '+') {
    qStart = qPosList[i];
    qEnd = qPosList[j];
  }
  else {
    qSize = hashIntVal(qSizes, axt->qName);
    qStart = qSize - qPosList[j] + 1;
    qEnd = qSize - qPosList[i] + 1;
  }

  /* output */
  fprintf(tr->outFile, "%s\t%d\t%d\t%s\t%d\t%d\t%c\t%.2f\t", 
	  axt->tName, tPosList[i]-1, tPosList[j],
	  axt->qName, qStart-1, qEnd,
	  axt->qStrand, 100.0 * score / (j-i+1));
  printCigarString(tr->outFile, axt, i, j);
  fputs("\n", tr->outFile);
}


void scanAxt(struct axt *axt, struct hash *qSizes, struct hash *tFilterAll, struct hash *qFilterAll, struct slThreshold *thresholds)
/* Scan one axt alignment and print conserved elements found to stdout.
 * THIS IS THE CORE FUNCTION OF THIS PROGRAM.
 * Arguments:
 * axt - alignment
 * qSizes - query assembly chromosome sizes
 * tFilterAll, qFilterAll - index of regions to exclude from scan
 * winSize - size of sliding window
 * thresholds - linked list of thresholds to call CEs at, and corresponding output filehandles
 */
{
  /* Variables to keep track of things as we loop through the alignment */
  int i = 0; /* column in alignment */
  int tPos = axt->tStart; /* position in target sequence */
  int qPos = axt->qStart; /* position in query sequence */
  int nrColumns;  /* counter for nr of columns seen after mask */
  int score;    /* sliding window score */
  struct slThreshold *tr;

  /* Three arrays where each element corresponds to a column in the alignment. */
  int *profile = needLargeMem(axt->symCount * sizeof(int)); /* cumulative conservation profile */
  int *tPosList = needLargeMem(axt->symCount * sizeof(int)); /* target seq position */
  int *qPosList = needLargeMem(axt->symCount * sizeof(int)); /* query seq position */
  /* E.g. at alignment column 5, target position tPosList[4] is aligned with query position qPosList[4],
   *      and alignment columns 1-5 contain a total of profile[4] matches.
   * Note:
   *  - in these 3 arrays, elements that that correspond to masked regions are not set
   *  - profile[] begins from zero again after each mask
   *  - target and query positions are set to -1 at gaps. */
  
  /* tFilter and qFilter are pointers to sorted arrays of coordinate ranges that should not be scanned.
   * The calls to searchFilter find the first filter overlapping or following the alignment */
  struct range *tFilter = tFilterAll ? searchFilter(tFilterAll, axt->tName, axt->tStart+1) : NULL;
  struct range *qFilter = qFilterAll ? searchFilter(qFilterAll, axt->qName, axt->qStart+1) : NULL;

  /* Initialize CE bounds for each threshold */
  for(tr = thresholds; tr != NULL; tr = tr->next) {
    tr->ceStart = -1; /* set to -1 = no CE found */
  }

  /* Main loop: go through alignment */ 
  while(i < axt->symCount) { /* loop until we have looked at entire alignment */

    /* if inside a mask, fast forward past it */
    do {
      if(tFilter) {
	while(tFilter->end <= tPos) tFilter++;
	if(tFilter->start <= tPos) {
	  if(tFilter->end >= axt->tEnd) goto endScan; /* using goto to break out of nested loop */
	  while(tFilter->end > tPos) {
	    if(axt->tSym[i] != '-') tPos++;
	    if(axt->qSym[i] != '-') qPos++;
	    i++;
	  }
	  tFilter++;
	}
      }
      if(qFilter) {
	while(qFilter->end <= qPos) qFilter++;
	if(qFilter->start <= qPos) {
	  if(qFilter->end >= axt->qEnd) goto endScan; /* using goto to break out of nested loop */
	  while(qFilter->end > qPos) {
	    if(axt->tSym[i] != '-') tPos++;
	    if(axt->qSym[i] != '-') qPos++;
	    i++;
	  }
	  qFilter++;
	}
      }
    } while(tFilter && tFilter->start <= tPos);

    /* handle first position after mask */
    profile[i] = bpScores[ (int) axt->qSym[i] ][ (int) axt->tSym[i] ];
    tPosList[i] = axt->tSym[i] == '-' ? -1 : ++tPos;
    qPosList[i] = axt->qSym[i] == '-' ? -1 : ++qPos;
    nrColumns = 1;

    /* handle remaining positions */
    for(i++; i < axt->symCount; i++) {
      /* break out of loop if we have come to a mask */
      if((tFilter && tFilter->start <= tPos) || (qFilter && qFilter->start <= qPos)) break;
      /* set positions */
      tPosList[i] = axt->tSym[i] == '-' ? -1 : ++tPos;
      qPosList[i] = axt->qSym[i] == '-' ? -1 : ++qPos;
      /* set profile */
      profile[i] = profile[i-1] + bpScores[ (int) axt->qSym[i]][ (int) axt->tSym[i] ];
      /* increment nr of columns seen after mask */
      nrColumns++;
      /* loop over user-defined thresholds */
      for(tr = thresholds; tr != NULL; tr = tr->next) {
	/* if have have seen enough columns to cover a window, evaluate that window */
	if(nrColumns >= tr->winSize) {
	  /* compute and check window score */
	  score = nrColumns > tr->winSize ? profile[i] - profile[i - tr->winSize] : profile[i];
	  if(score >= tr->minScore) {
	    /* score is above threshold: begin or extend conserved element */
	    if(tr->ceStart == -1) tr->ceStart = i - tr->winSize + 1;
	    tr->ceEnd = i;
	  }
	  else {
	    /* score is below threshold: close and print any open conserved elements that are more than a window away */
	    if(tr->ceStart != -1 && tr->ceEnd < i - tr->winSize + 1) {
	      printElement(tr, axt, qSizes, profile, tPosList, qPosList);
	      tr->ceStart = -1;
	    }
	  }
	}  
      }
    }

    /* close and print any open conserved elements */
    for(tr = thresholds; tr != NULL; tr = tr->next) {
      if(tr->ceStart != -1) {
	printElement(tr, axt, qSizes, profile, tPosList, qPosList);
	tr->ceStart = -1;
      }
    }
  }
  
 endScan:

  /* free memory */
  freez(&profile);
  freez(&tPosList);
  freez(&qPosList);
}


void ceScan(char *tFilterFile, char *qFilterFile, char *qSizeFile, struct slThreshold *thresholds, int nrAxtFiles, char *axtFiles[])
/* ceScan - Find conserved elements. */
{
  struct lineFile *lf;
  struct axt *axt;
  struct hash *tFilter, *qFilter, *qFilterRev, *qSizes;
  int i;

  setBpScores(bpScores);
  qSizes = loadIntHash(qSizeFile);
  tFilter = tFilterFile ? readFilter(tFilterFile) : NULL;
  qFilter = qFilterFile ? readFilter(qFilterFile) : NULL;
  qFilterRev = qFilter ? makeReversedFilter(qFilter, qSizes) : NULL;

  for(i = 0; i < nrAxtFiles; i++) {
    lf = lineFileOpen(axtFiles[i], TRUE);
    while ((axt = axtRead(lf)) != NULL) {     
      scanAxt(axt,
	      qSizes,
	      tFilter,
	      axt->qStrand == '+' ? qFilter : qFilterRev,
	      thresholds);
      axtFree(&axt);
    }
    lineFileClose(&lf);
  }

  /* here we should free the hashes allocated above if this function was part of a larger program */
}


int main(int argc, char *argv[])
/* Process command line. */
{
  char *tFilterFile, *qFilterFile, *qSizeFile, *outFilePrefix;
  int i, n, minScore, winSize;
  struct slThreshold *trList = NULL, *tr;
  char rest, path[PATH_LEN];

  /* Parse options */
  optionInit(&argc, argv, options);
  tFilterFile = optionVal("tFilter", NULL);
  qFilterFile = optionVal("qFilter", NULL);
  outFilePrefix = optionVal("outPrefix", NULL);

  /* Parse mandatory arguments */
  if (argc < 4) usage();

  /* First arg is a string */
  qSizeFile = argv[1];

  /* Argument 4 and possibly some more are threshold specifiers. */
  /* Parse threshold specifiers and open output files. */
  for(i = 2; i < argc-1; i++) {
    n = sscanf(argv[i],"%d,%d%c", &minScore, &winSize, &rest);
    if(n != 2) {
      if(trList == NULL) errAbort("Invalid threshold string: %s", argv[i]);
      break;
    }
    tr = needMem(sizeof(*tr));
    tr->minScore = minScore;
    tr->winSize = winSize;   
    if(outFilePrefix) {
      safef(path, sizeof(path), "%s_%d_%d", outFilePrefix, minScore, winSize);
      tr->outFile = mustOpen(path, "w");
    }
    else tr->outFile = stdout;
    slAddHead(&trList, tr);
  }

  /* The remanining arguments (starting at argv+i) should be names of axt files */

  /* Call function ceScan with the arguments */
  ceScan(tFilterFile, qFilterFile, qSizeFile, trList, argc-i, argv+i);

  /* Close all output files */
  if(outFilePrefix)
    for(tr = trList; tr != NULL; tr = tr->next)
      fclose(tr->outFile);

  /* Return success */
  return 0;
}
