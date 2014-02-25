#include "CNEr.h"


/*############################################*/

SEXP myReadBed(SEXP filepath){
  // load a filter file into R, and to be a GRanges in 1-based coordinates
  // This is tested and without memory leak!
  filepath = AS_CHARACTER(filepath);
  if(!IS_CHARACTER(filepath) || LENGTH(filepath) != 1)
    error("'filepath' must be a single string");
  if(STRING_ELT(filepath, 0) == NA_STRING)
    error("'filepath' is NA");
  // If filepath_elt is defined this way, the memory will be reclaimed by the end of .Call by R, do not need the free()
  char *filepath_elt = R_alloc(strlen(CHAR(STRING_ELT(filepath, 0))) + 1, sizeof(char));
  strcpy(filepath_elt, CHAR(STRING_ELT(filepath, 0)));
  Rprintf("Reading %s \n", filepath_elt);
  struct lineFile *lf = lineFileOpen(filepath_elt, TRUE);
  char *row[3];
  int nRanges = 0;
  while(lineFileRow(lf, row)){
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    nRanges++;
  }
  lineFileClose(&lf);
  SEXP chromNames, starts, ends, returnList;
  PROTECT(returnList = NEW_LIST(3));
  chromNames = NEW_CHARACTER(nRanges);
  SET_VECTOR_ELT(returnList, 0, chromNames);
  starts = NEW_INTEGER(nRanges);
  SET_VECTOR_ELT(returnList, 1, starts);
  ends = NEW_INTEGER(nRanges);
  SET_VECTOR_ELT(returnList, 2, ends);
  int *p_starts, *p_ends;
  int j = 0;
  p_starts = INTEGER_POINTER(starts);
  p_ends = INTEGER_POINTER(ends);
  lf = lineFileOpen(filepath_elt, TRUE);
  while(lineFileRow(lf, row)){
    if(sameString(row[0], "track") || sameString(row[0], "browser")) continue;
    p_starts[j] = lineFileNeedNum(lf, row, 1) + 1;
    p_ends[j] = lineFileNeedNum(lf, row, 2);
    if(p_starts[j] > p_ends[j])
      errAbort("start after end line %d of %s", lf->lineIx, lf->fileName);
    SET_STRING_ELT(chromNames, j, mkChar(row[0]));
    j++;
  }
  lineFileClose(&lf);
  UNPROTECT(1);
  return(returnList);
}


/* ------------------------- Use the DNAStringSet way--------------------------*/
// The advantage of this way is that it allocates the large memory together (in Linux, the mmep way) which can be released by the OS immediately after deletion.
SEXP axt_info(SEXP filepath){
  // read a axt file and get the alignment length
  filepath = AS_CHARACTER(filepath);
  int nrAxtFiles, i;
  nrAxtFiles = GET_LENGTH(filepath);
  Rprintf("The number of axt files %d\n", nrAxtFiles);
  struct axt *curAxt;
  struct lineFile *lf;
  IntAE width_buf;
  width_buf = new_IntAE(0, 0, 0);
  char *filepath_elt;
  for(i = 0; i < nrAxtFiles; i++){
    filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, i)))+1, sizeof(char));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, i)));
    lf = lineFileOpen(filepath_elt, TRUE);
    while((curAxt = axtRead(lf)) != NULL){
      IntAE_insert_at(&width_buf, IntAE_get_nelt(&width_buf), curAxt->symCount);
      axtFree(&curAxt);
    }
    lineFileClose(&lf);
  }
  SEXP width;
  PROTECT(width = new_INTEGER_from_IntAE(&width_buf));
  Rprintf("The number of axt alignments is %d\n", GET_LENGTH(width));
  UNPROTECT(1);
  return(width);
}

SEXP readAxt(SEXP filepath){
  // load a axt file into R, and to be axt object
  // This is tested and without memory leak!
  // The jim kent's axt struct holds the starts in 0-based. Here we put it into 1-based.
  filepath = AS_CHARACTER(filepath);
  int nrAxtFiles, i, nrAxts;
  nrAxtFiles = GET_LENGTH(filepath);
  struct axt *axt=NULL, *curAxt;
  struct lineFile *lf;
  IntAE width_buf;
  SEXP ans_qSym, ans_tSym, width;
  SEXP qNames, qStart, qEnd, qStrand, tNames, tStart, tEnd, tStrand, score, symCount, returnList;
  PROTECT(returnList = NEW_LIST(12));
  PROTECT(width = axt_info(filepath));
  nrAxts = GET_LENGTH(width);
  tNames = NEW_CHARACTER(nrAxts);
  SET_VECTOR_ELT(returnList, 0, tNames);
  tStart = NEW_INTEGER(nrAxts);
  SET_VECTOR_ELT(returnList, 1, tStart);
  tEnd = NEW_INTEGER(nrAxts);
  SET_VECTOR_ELT(returnList, 2, tEnd);
  tStrand = NEW_CHARACTER(nrAxts);
  SET_VECTOR_ELT(returnList, 3, tStrand);
  ans_tSym = alloc_XRawList("BStringSet", "BString", width);
  SET_VECTOR_ELT(returnList, 4, ans_tSym);
  qNames = NEW_CHARACTER(nrAxts);
  SET_VECTOR_ELT(returnList, 5, qNames);
  qStart = NEW_INTEGER(nrAxts);
  SET_VECTOR_ELT(returnList, 6, qStart);
  qEnd = NEW_INTEGER(nrAxts);
  SET_VECTOR_ELT(returnList, 7, qEnd);
  qStrand = NEW_CHARACTER(nrAxts);
  SET_VECTOR_ELT(returnList, 8, qStrand);
  ans_qSym = alloc_XRawList("BStringSet", "BString", width);
  SET_VECTOR_ELT(returnList, 9, ans_qSym);
  score = NEW_INTEGER(nrAxts);
  SET_VECTOR_ELT(returnList, 10, score);
  symCount = NEW_INTEGER(nrAxts);
  SET_VECTOR_ELT(returnList, 11, symCount);
  UNPROTECT(1);
  cachedXVectorList cached_ans_qSym, cached_ans_tSym;
  cachedCharSeq cached_ans_elt;
  cached_ans_qSym = cache_XVectorList(ans_qSym);
  cached_ans_tSym = cache_XVectorList(ans_tSym);
  int *p_qStart, *p_qEnd, *p_tStart, *p_tEnd, *p_score, *p_symCount;
  p_qStart = INTEGER_POINTER(qStart);
  p_qEnd = INTEGER_POINTER(qEnd);
  p_tStart = INTEGER_POINTER(tStart);
  p_tEnd = INTEGER_POINTER(tEnd);
  p_score = INTEGER_POINTER(score);
  p_symCount = INTEGER_POINTER(symCount);
  int j = 0;
  i = 0;
  for(j = 0; j < nrAxtFiles; j++){
    char *filepath_elt = (char *) R_alloc(strlen(CHAR(STRING_ELT(filepath, j)))+1, sizeof(char));
    strcpy(filepath_elt, CHAR(STRING_ELT(filepath, j)));
    lf = lineFileOpen(filepath_elt, TRUE);
    while((axt = axtRead(lf)) != NULL){
      SET_STRING_ELT(qNames, i, mkChar(axt->qName));
      p_qStart[i] = axt->qStart + 1;
      p_qEnd[i] = axt->qEnd;
      if(axt->qStrand == '+')
        SET_STRING_ELT(qStrand, i, mkChar("+"));
      else
        SET_STRING_ELT(qStrand, i, mkChar("-"));
      cached_ans_elt = get_cachedXRawList_elt(&cached_ans_qSym, i);
      memcpy((char *) (&cached_ans_elt)->seq, axt->qSym, axt->symCount * sizeof(char));

      SET_STRING_ELT(tNames, i, mkChar(axt->tName));
      p_tStart[i] = axt->tStart + 1;
      p_tEnd[i] = axt->tEnd;
      if(axt->tStrand == '+')
        SET_STRING_ELT(tStrand, i, mkChar("+"));
      else
        SET_STRING_ELT(tStrand, i, mkChar("-"));
      cached_ans_elt = get_cachedXRawList_elt(&cached_ans_tSym, i);
      memcpy((char *) (&cached_ans_elt)->seq, axt->tSym, axt->symCount * sizeof(char));
      p_score[i] = axt->score;
      p_symCount[i] = axt->symCount;
      i++;
      axtFree(&axt);
    }
    lineFileClose(&lf);
  }
  UNPROTECT(1);
  //return R_NilValue;
  return returnList;
}

