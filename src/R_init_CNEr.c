#include "CNEr.h"

#define CALLMETHOD_DEF(fun, numArgs) {#fun, (DL_FUNC) &fun, numArgs}

static const R_CallMethodDef callMethods[] = {
  /* io.c */
  CALLMETHOD_DEF(myReadBed, 1),
  CALLMETHOD_DEF(myReadAxt, 1),
  CALLMETHOD_DEF(axt_info, 1),
  
  /* utils.c */
  CALLMETHOD_DEF(bin_from_coord_range, 2),
  CALLMETHOD_DEF(bin_ranges_from_coord_range, 2),
  
  /* ceScan.c */
  CALLMETHOD_DEF(myCeScan, 23),
  CALLMETHOD_DEF(myCeScanFile, 8),
  {NULL, NULL, 0}
};

void R_init_CNEr(DllInfo *info)
{
  R_registerRoutines(info, NULL, callMethods, NULL, NULL);
  return;
}

