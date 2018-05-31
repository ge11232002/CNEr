/********************************************************************
 *                        Alignment related function
 *                        Author: Ge Tan
 *******************************************************************/
#include "Rdefines.h"

SEXP subAlignment(SEXP axttStart, SEXP axttEnd, SEXP axttSym,
                  SEXP axtqStart, SEXP axtqEnd, SEXP axtqSym,
                  SEXP subtStart, SEXP subtEnd, SEXP symCount){
// The input must have same length and subtStart,subtEnd must be inside 
// the range of axttStart,axttEnd.
  int nrAxt = GET_LENGTH(symCount);
  int i, j;
  SEXP returnList;
  PROTECT(returnList = NEW_LIST(6));
  SEXP newtStart, newtEnd, newtSym, newqStart, newqEnd, newqSym;
  newtStart = NEW_INTEGER(nrAxt);
  newtEnd = NEW_INTEGER(nrAxt);
  newtSym = NEW_CHARACTER(nrAxt);
  newqStart = NEW_INTEGER(nrAxt);
  newqEnd = NEW_INTEGER(nrAxt);
  newqSym = NEW_CHARACTER(nrAxt);
  SET_VECTOR_ELT(returnList, 0, newtStart);
  SET_VECTOR_ELT(returnList, 1, newtEnd);
  SET_VECTOR_ELT(returnList, 2, newtSym);
  SET_VECTOR_ELT(returnList, 3, newqStart);
  SET_VECTOR_ELT(returnList, 4, newqEnd);
  SET_VECTOR_ELT(returnList, 5, newqSym);
  int nrGapsTarget, nrGapsQuery, cpStartTarget, cpEndTarget, cpStartQuery, cpEndQuery;
  for(i=0; i<nrAxt; i++){
    // Rprintf("The char %d \n", i);
    nrGapsTarget = 0;
    nrGapsQuery = 0;
    cpStartTarget = 0;
    cpEndTarget = 0;
    for(j=0; j<INTEGER(symCount)[i]; j++){
      if(CHAR(STRING_ELT(axttSym, i))[j] == '_' || 
          CHAR(STRING_ELT(axttSym, i))[j] == '-'){
        //Rprintf("%c ", CHAR(STRING_ELT(axttSym, i))[j]);
        nrGapsTarget++;
      }
      if(CHAR(STRING_ELT(axtqSym, i))[j] == '_' ||
          CHAR(STRING_ELT(axtqSym, i))[j] == '-'){
        nrGapsQuery++;
      }
      if(INTEGER(subtStart) == INTEGER(axttStart) + j - nrGapsTarget){
        cpStartTarget = j;
        INTEGER(newqStart) = j + INTEGER(axtqStart);
      }
      if(INTEGER(subtEnd) == INTEGER(axttStart) + j - nrGapsTarget){
        cpEndTarget = j;
      }
    }
    Rprintf("\n");
  }
  UNPROTECT(1);
  return(R_NilValue);
}




