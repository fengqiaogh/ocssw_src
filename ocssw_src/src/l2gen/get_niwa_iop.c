/*
 * get_niwa_iop.c -- wrapper functions to interface NIWA/UoP/Moore IOP algorithm
 *                   to prodgen and convl12
 */

#include "l12_proto.h"
#include "l2prod.h"

#include "niwa_iop.h"


/* keep track of last scan line processed */
static int last_recnum = -1;


/* temporary arrays for iop values for this scan line */
static float *niwa_a; /* absorption coefficient per band and pixel */
static float *niwa_bb; /* backscatter coefficient per band and pixel */
static int16 *niwa_iopf; /* iop flags per pixel */

/*
 * check if we have already calculated the values for this scan line
 */
static int niwa_ran(int recnum) {
    return recnum == last_recnum;
}

/* 
 * allocate space for one scan line of IOP variables 
 */
static void alloc_niwa(int32_t npix, int32_t nbands) {
    if ((niwa_a = calloc(npix * nbands, sizeof (float))) == NULL) {
        fprintf(stderr, "-E- %s line %d : error allocating memory for NIWA IOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((niwa_bb = calloc(npix * nbands, sizeof (float))) == NULL) {
        fprintf(stderr, "-E- %s line %d : error allocating memory for NIWA IOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if ((niwa_iopf = calloc(npix, sizeof (int16))) == NULL) {
        fprintf(stderr, "-E- %s line %d : error allocating memory for NIWA IOP.\n",
                __FILE__, __LINE__);
        exit(1);
    }
}

static void run_niwa(l2str *l2rec) {
    static int first_time = 1;

    if (first_time) {
        /* allocate a, bb & flag arrays */
        alloc_niwa(l2rec->l1rec->npix, l2rec->l1rec->l1file->nbands);
        first_time = 0;
    }

    /* call the NIWA iop code to calculate the a and bb values */
    niwa_iop(l2rec, niwa_a, niwa_bb, niwa_iopf);

    last_recnum = l2rec->l1rec->iscan;
}

/*
 * interface to prodgen() for NIWA iops
 */
void get_niwa(l2str *l2rec, l2prodstr *p, float prod[]) {
    int prod_id = p->cat_ix;
    int ib = p->prod_ix;
    int ip, ipb;

    if (!niwa_ran(l2rec->l1rec->iscan))
        run_niwa(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        ipb = ip * l2rec->l1rec->l1file->nbands + ib;

        switch (prod_id) {
        case CAT_a_niwa:
            prod[ip] = niwa_a[ipb];
            break;

        case CAT_bb_niwa:
            prod[ip] = niwa_bb[ipb];
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to NIWA IOP\n",
                    __FILE__, __LINE__, prod_id);
            exit(1);
        }
    }
}

/*
 * get the NIWA iop status flags 
 */
int16 * get_flags_niwa(l2str *l2rec) {
    if (!niwa_ran(l2rec->l1rec->iscan))
        run_niwa(l2rec);

    return niwa_iopf;
}

/* 
 * interface to convl12() for NIWA iops
 */
void iops_niwa(l2str *l2rec) {
    int ib, ip, ipb;

    l1str *l1rec = l2rec->l1rec;
    int32_t nbands = l1rec->l1file->nbands;

    if (!niwa_ran(l1rec->iscan))
        run_niwa(l2rec);

    for (ip = 0; ip < l1rec->npix; ip++) {
        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            l2rec->a[ipb] = niwa_a[ipb];
            l2rec->bb[ipb] = niwa_bb[ipb];
        }
    }
}

