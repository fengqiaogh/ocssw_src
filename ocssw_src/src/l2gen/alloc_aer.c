#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "aer_struc.h"
#include "l12_parms.h"

/* --------------------------------------------------------- */
/* alloc_aer() - allocates 1 aerosol record to hold data     */
/*              for a single scan of "npix" pixels.          */

/* --------------------------------------------------------- */
int alloc_aer(int32_t npix, int32_t nbands, aestr *rec) {
    char *p;
    int32_t len = 2 * sizeof (int32_t) * npix /* model #      */
            + 1 * sizeof (float)*npix /* model ratio  */
            + 1 * nbands * sizeof (float)*npix; /* taua         */

    if ((p = (char *) malloc(len)) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Memory allocation failure.\n",
                __FILE__, __LINE__);
        return (0);
    }
    memset(p, '\0', len);

    rec->length = len;
    rec->npix = npix;
    rec->data = p;

    rec->mod_min = (int32_t *) p;
    p += sizeof (int32_t) * npix;
    rec->mod_max = (int32_t *) p;
    p += sizeof (int32_t) * npix;
    rec->mod_rat = (float *) p;
    p += sizeof (float)*npix;
    rec->taua = (float *) p;
    p += sizeof (float)*npix*nbands;

    return (len);
}
