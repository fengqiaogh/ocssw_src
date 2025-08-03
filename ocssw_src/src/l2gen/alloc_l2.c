#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include "l2_struc.h"
#include "l12_proto.h"

void free_l2(l2str *l2rec) {
    free((void *) l2rec->data);
}


/* --------------------------------------------------------- */
/* alloc_l2() - allocates 1 level-2 record to hold data for  */
/*              a single scan of "npix" pixels.              */

/* --------------------------------------------------------- */
int alloc_l2(l1str *l1rec, l2str *l2rec) {
    l2rec->l1rec = l1rec;

    int32_t npix = l1rec->npix;
    int32_t nbands = l1rec->l1file->nbands;
    char *p;

    if(input->proc_uncertainty){
        if(input->aer_opt!=AERRHMSEPS && input->aer_opt!=AERRHSM){
            printf("-E- %s line %d: Unable to calculate the uncertainty for aer_opt other than -17 or -18.\n", __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }

        l1rec->uncertainty= (uncertainty_t*) malloc(sizeof (uncertainty_t));
        if (alloc_uncertainty(nbands, input->nbands_ac, npix, l1rec->uncertainty) != 0) {
            printf("-E- %s line %d: Unable to allocate error record.\n", __FILE__, __LINE__);
            exit(FATAL_ERROR);
        }
    }

    int32_t len = 5 * sizeof (int32_t) * npix
            + 7 * sizeof (float)*npix
            + 10 * sizeof (float)*npix * nbands;

    if(l1rec->uncertainty){
        len+=sizeof (float)*npix * nbands;   // for  Rrs_unc
        if(input->proc_uncertainty==2)
            len+=sizeof(float)*npix*nbands*nbands;   //for covariance matrix
    }


    if (len % 4 != 0)
        len = len / 4 * 4 + 4;

    if ((p = (char *) malloc(len)) == NULL) {
        fprintf(stderr,
                "-E- %s line %d: Memory allocation failure.\n",
                __FILE__, __LINE__);
        return (0);
    }

    l2rec->length = len;
    l2rec->data = p;

    // var[npix]
    l2rec->num_iter = (int32_t *) p;
    p += sizeof(int32_t)*npix;
    l2rec->aermodmin = (int32_t *) p;
    p += sizeof(int32_t)*npix;
    l2rec->aermodmax = (int32_t *) p;
    p += sizeof(int32_t)*npix;
    l2rec->aermodmin2 = (int32_t *) p;
    p += sizeof(int32_t)*npix;
    l2rec->aermodmax2 = (int32_t *) p;
    p += sizeof(int32_t)*npix;

    l2rec->chl = (float *) p;
    p += sizeof(float)*npix;
    l2rec->eps = (float *) p;
    p += sizeof(float)*npix;
    l2rec->chi2 = (float *) p;
    p += sizeof(float)*npix;
    l2rec->aerratio = (float *) p;
    p += sizeof(float)*npix;
    l2rec->aerratio2 = (float *) p;
    p += sizeof(float)*npix;
    l2rec->aerindex = (float *) p;
    p += sizeof(float)*npix;

    // var[npix][nbands]
    l2rec->taua = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->La = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->Lw = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->nLw = (float *) p;
    p += sizeof(float)*npix*nbands;

    l2rec->brdf = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->Rrs = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->Rrs_raman = (float *) p;
    p += sizeof(float)*npix*nbands;
    if(l1rec->uncertainty){
        l2rec->Rrs_unc = (float *) p;
        p += sizeof (float)*npix*nbands;
    } else {
        l2rec->Rrs_unc = NULL;
    }
    l2rec->outband_correction = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->a = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->bb = (float *) p;
    p += sizeof(float)*npix*nbands;
    l2rec->chl_unc = (float *) p;
    p += sizeof (float)*npix;

    if(input->proc_uncertainty==2){
        l2rec->covariance_matrix = (float *) p;
        p += sizeof (float) * npix*nbands*nbands;
    }else {
        l2rec->covariance_matrix = NULL;
    }

    if ((len - (int32_t) (p - l2rec->data)) < 0) {
        printf("%s Line %d: bad allocation on L2 record\n", __FILE__, __LINE__);
        exit(1);
    }

    // init to NULL
    l2rec->bindx = NULL;
    l2rec->sst = NULL;
    l2rec->tgrec = NULL;

    printf("Allocated %d bytes in L2 record.\n", (int) (p - l2rec->data));

    return (len);
}
