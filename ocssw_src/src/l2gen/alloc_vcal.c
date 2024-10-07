#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "vcal_struc.h"
#include "l12_parms.h"

/* --------------------------------------------------------- */
/* alloc_target() - allocates 1 target record to hold data   */
/*              for a single scan of "npix" pixels.          */

/* --------------------------------------------------------- */
int alloc_vcal(int32_t npix, int32_t nbands, vcstr *rec) {
    char *p;
    int32_t len = sizeof (float)*npix * nbands * 6;

    if ((p = (char *) malloc(len)) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Memory allocation failure.\n",
                __FILE__, __LINE__);
        return (0);
    }
    memset(p, '\0', len);

    rec->length = len;
    rec->npix = npix;
    rec->nbands = nbands;
    rec->data = p;

    rec->vLt = (float *) p;
    p += sizeof (float)*npix*nbands;
    rec->Lw = (float *) p;
    p += sizeof (float)*npix*nbands;
    rec->nLw = (float *) p;
    p += sizeof (float)*npix*nbands;
    rec->tLw = (float *) p;
    p += sizeof (float)*npix*nbands;
    rec->brdfsat = (float *) p;
    p += sizeof (float)*npix*nbands;
    rec->brdftgt = (float *) p;
    p += sizeof (float)*npix*nbands;

    return (len);
}
