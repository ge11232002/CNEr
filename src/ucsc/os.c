/* Added by rtracklayer to condition on the platform */
/* License: Artistic-2.0 */

#ifdef WIN32
#include "oswin9x.c"
#else
#include "osunix.c"
#endif
