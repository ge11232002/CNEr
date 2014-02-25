#include "CNEr.h"

// This is borrowed from Jim Kend's source, binRange.c.
/* add one new level to get coverage past chrom sizes of 512 Mb
 *  *  effective limit is now the size of an integer since chrom start
 *   *  and end coordinates are always being used in int's == 2Gb-1 */
static int binOffsetsExtended[] = 
  {4096+512+64+8+1, 512+64+8+1, 64+8+1, 8+1, 1, 0};

static int binOffsets[] = 
  {512+64+8+1, 64+8+1, 8+1, 1, 0};

SEXP bin_from_coord_range(SEXP starts, SEXP ends){
//Return the bin numbers that should be assigned to a feature 
//spanning the given ranges.
//Here the inputs of starts and ends are 1-based
  starts = AS_INTEGER(starts);
  ends = AS_INTEGER(ends);
  int i, n;
  int *p_start, *p_end, *p_bin;
  n = GET_LENGTH(starts);
  SEXP bins;
  PROTECT(bins = NEW_INTEGER(n));
  p_start = INTEGER_POINTER(starts);
  p_end = INTEGER_POINTER(ends);
  p_bin = INTEGER_POINTER(bins);
  for(i=0; i<n; i++){
    // passing the coordinates to binFromRange in 0-based for start, 
    // 1-based for end.
    p_bin[i] = binFromRange(p_start[i]-1, p_end[i]);
  }
  UNPROTECT(1);
  return bins;
}


SEXP bin_ranges_from_coord_range_standard(SEXP start, SEXP end){
  // both are 1-based
  if(GET_LENGTH(start) != 1 || GET_LENGTH(end) != 1)
    error("'start' and 'end' must be a single integer");
  start = AS_INTEGER(start);
  end = AS_INTEGER(end);
  int startBin = INTEGER(start)[0] - 1;
  int endBin = INTEGER(end)[0] - 1;
  int _binFirstShift = binFirstShift();
  int _binNextShift = binNextShift();
  startBin >>= binFirstShift();
  endBin >>= binFirstShift();
  int n = ArraySize(binOffsets);
  SEXP ans;
  PROTECT(ans = allocMatrix(INTSXP, n, 2));
  int *rans = INTEGER_POINTER(ans);
  // we use a matrix to store the bin ranges for each bin level.
  int i;
  for(i=0; i<n; ++i){
    rans[i] = binOffsets[i] + startBin;
    rans[i+n] = binOffsets[i] + endBin;
    startBin >>= binNextShift();
    endBin >>= binNextShift();
  }
  UNPROTECT(1);
  return ans;
}

SEXP bin_ranges_from_coord_range_extended(SEXP start, SEXP end){
  // both are 1-based
  if(GET_LENGTH(start) != 1 || GET_LENGTH(end) != 1)
    error("'start' and 'end' must be a single integer");
  start = AS_INTEGER(start);
  end = AS_INTEGER(end);
  int startBin = INTEGER(start)[0] - 1;
  int endBin = INTEGER(end)[0] - 1;
  startBin >>= binFirstShift();
  endBin >>= binFirstShift();
  int n = ArraySize(binOffsetsExtended);
  SEXP ans;
  PROTECT(ans = allocMatrix(INTSXP, n, 2));
  int *rans = INTEGER_POINTER(ans);
  int i;
  for(i=0; i<n; ++i){
    rans[i] = _binOffsetOldToExtended + binOffsetsExtended[i] + startBin;
    rans[i+n] = _binOffsetOldToExtended + binOffsetsExtended[i] + endBin;
    startBin >>= binNextShift();
    endBin >>= binNextShift();
  }
  UNPROTECT(1);
  return ans;
}

SEXP bin_ranges_from_coord_range(SEXP start, SEXP end){
  //Return the set of bin ranges that overlap a given coordinate range. 
  //It is usually more convenient to use bin_restriction string than 
  //to use this method directly.
  // Here, start and end are 1-based.
  if (INTEGER(AS_INTEGER(end))[0] <= BINRANGE_MAXEND_512M)
    return bin_ranges_from_coord_range_standard(start, end);
  else
    return bin_ranges_from_coord_range_extended(start, end);
}



/*SEXP calc_window_scores(SEXP coverage, SEXP context_start, SEXP context_end, SEXP win_nr_steps, SEXP step_size){
  // Here the start and end are 1-based
  win_nr_steps = AS_INTEGER(win_nr_steps);
  step_size = AS_INTEGER(step_size);
  context_start = AS_INTEGER(context_start);
  context_end = AS_INTEGER(context_end);
  int context_size = INTEGER(context_end)[0] - INTEGER(context_start)[0] + 1;
  int nr_blocks;
  nr_blocks = (((context_size % INTEGER(step_size)[0]) != 0)?(context_size/INTEGER(step_size)[0] + 1):(context_size/INTEGER(step_size)[0]));
  Rprintf("The nr_blocks %d\n", nr_blocks);
  //int blk_scores[((nr_blocks>INTEGER(win_nr_steps)[0])?(nr_blocks):(INTEGER(win_nr_steps)[0]+1))] = {0};
  int *blk_scores;
  blk_scores = (int *) R_alloc(((nr_blocks>INTEGER(win_nr_steps)[0])?(nr_blocks):(INTEGER(win_nr_steps)[0]+1)), sizeof(int));
  int blk_start = INTEGER(context_start)[0];
  int blk_end = INTEGER(context_start)[0] + INTEGER(step_size)[0] - 1;
  int prev_win_score = 0; 
  int i;// The blocks index
  for(i=0; i < nr_blocks

  return R_NilValue;
}*/
