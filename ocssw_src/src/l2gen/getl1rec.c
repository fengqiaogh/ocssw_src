#include <stdio.h>
#include <stdlib.h>
#include "l1.h"
#include "l12_parms.h"
#include "l12_proto.h"

l1qstr l1que;

/* --------------------------------------------------------- */
/* free_l1q() - free memory allocated to the l1 queue        */

/* --------------------------------------------------------- */
void free_l1q() {
    int32_t i;
    int32_t nq = l1que.nq;

    if (nq > 0)
        for (i = 0; i < nq; i++)
            free_l1(&(l1que.r[i]));
}


/* --------------------------------------------------------- */
/* alloc_l1q() - allocates a structure to hold nq level-1b   */
/*   records                                                 */

/* --------------------------------------------------------- */
int32_t alloc_l1q(int32_t nq, filehandle *l1file, l1qstr *l1que) {
    int32_t i;

    if (nq <= 0)
        nq = 1;

    if (nq > NQMAX) {
        printf("Queue size limit of %d exceeded: %d\n", NQMAX, nq);
        return (0);
    }

    /* force the que size to an odd number and init the center scan num */
    l1que->nq = (nq / 2)*2 + 1;
    l1que->cscan = -1;

    for (i = 0; i < nq; i++) {
        if (alloc_l1(l1file, &(l1que->r[i])) == 0) {
            fprintf(stderr,
                    "-E- %s line %d: Memory allocation failure at L1B scan %d.\n",
                    __FILE__, __LINE__, i);
            return (0);
        }
    }


    return (nq);
}


/* ---------------------------------------------------------------- */
/* Read one or more level 1 records from the file pointed to by the */
/* input file handle, as required to produce a queue of l1 records  */
/* centered on the input scan number.                               */

/* ---------------------------------------------------------------- */
int loadl1q(filehandle *l1file, int32_t iscan, int32_t dscan) {
    int32_t nq = l1que.nq;
    int32_t i, iq, recnum;
    int32_t fscan = iscan - nq / 2;
    int32_t lscan = iscan + nq / 2;
    l1str tmp;

    if (iscan < 0 || iscan > l1file->nscan)
        return (EXIT_FAILURE);
    /*                                                              */
    /* If the current queue center is not the preceeding scan, then */
    /* this is the first call or we are reading out of sequence, so */
    /* we will just load the entire queue from the input file.      */
    /* Otherwise, we just need to shift the queue and read a scan.  */
    /*                                                              */
    if (l1que.cscan == -1 || l1que.cscan != iscan - dscan) {

        for (i = fscan; i <= lscan; i++) {
            iq = i - fscan;
            recnum = MIN(MAX(MAX(i, l1_input->sline - 1), 0), l1file->nscan - 1);
            if (readl1(l1file, recnum, &l1que.r[iq]) != 0)
                return (EXIT_FAILURE);
            if (loadl1(l1file, &l1que.r[iq]) != 0)
                return (EXIT_FAILURE);
        }

    } else {

        /*                                                          */
        /* We really just re-arrange the record data pointers       */
        /*                                                          */
        memcpy(&tmp, &l1que.r[0], sizeof (l1str));
        for (iq = 1; iq < nq; iq++)
            memcpy(&l1que.r[iq - 1], &l1que.r[iq], sizeof (l1str));
        memcpy(&l1que.r[nq - 1], &tmp, sizeof (l1str));

        /*                                                          */
        /* Now read the next scan into the top of the queue         */
        /*                                                          */
        recnum = MIN(MAX(MAX(iscan + nq / 2, l1_input->sline - 1), 0), l1file->nscan - 1);
        if (readl1(l1file, recnum, &l1que.r[nq - 1]) != 0)
            return (EXIT_FAILURE);
        if (loadl1(l1file, &l1que.r[nq - 1]) != 0)
            return (EXIT_FAILURE);
    }

    /* We now have a queue of nq L1B records, centered on iscan     */
    l1que.cscan = iscan;

    return (EXIT_SUCCESS);
}

void lt_agregat_ocis(l1str*l1rec) {
    static int firstcall = 1;

    int i, j, ipb, k;
    float *wave = l1rec->l1file->fwave;
    int npix = l1rec->npix;
    int nbands = l1rec->l1file->nbands;
    float templt = 0;

    int ib325 = windex(325, wave, nbands); //5nm
    int ib340 = windex(340, wave, nbands); //5nm
    int ib665 = windex(665, wave, nbands); //10nm
    int ib680 = windex(680, wave, nbands); //10nm
    int ib750 = windex(750, wave, nbands); //10nm
    int ib865 = windex(865, wave, nbands); //40nm

    static float *Lt_scan;

    if (firstcall) {
        Lt_scan = (float*) malloc(npix * nbands * sizeof (float));
        firstcall = 0;
    }
    memcpy(Lt_scan, l1rec->Lt, sizeof (float)*npix * nbands);


    for (i = 0; i < npix; i++) {
        ipb = i*nbands;
        for (j = 0; j < nbands; j++) {
            if (j == ib325 || j == ib340)
                continue;
            else if (j == ib865) {
                templt = 0;
                for (k = j - 4; k <= j + 4; k++)
                    templt += Lt_scan[ipb + k];
                l1rec->Lt[ipb + j] = templt / 9.0;
            } else if (j == ib665 || j == ib680 || j == ib750) {
                l1rec->Lt[ipb + j] = (Lt_scan[ipb + j - 1] + Lt_scan[ipb + j]) / 2.0;
            } else {
                if (j ==0)
                    l1rec->Lt[ipb + j] = (2 * Lt_scan[ipb + j] + Lt_scan[ipb + j + 1]) / 3.0;                    
                else if (j == nbands - 1)
                    l1rec->Lt[ipb + j] = (Lt_scan[ipb + j - 1] + 2 * Lt_scan[ipb + j]) / 3.0;
                else
                    l1rec->Lt[ipb + j] = (Lt_scan[ipb + j - 1] + Lt_scan[ipb + j] + Lt_scan[ipb + j + 1]) / 3.0;
            }
        }
    }

}


/* ---------------------------------------------------------------- */

/* ---------------------------------------------------------------- */
int getl1rec(int32_t iscan, int32_t dscan, l1str *l1rec) {
    static int32_t npix = -1;
    static int32_t nq = -1;
    l1str *crec;
    int32_t ip,iq;
    filehandle* l1file = l1rec->l1file;
    int nbands=l1file->nbands;

    if (nq == -1)
        nq = MAX(input->fctl.nscan, NQMIN);

    crec = &l1que.r[nq / 2];

    /*                                                              */
    /* if the queue size is inconsistent with the filter scan count */
    /* then this must be the first time through, so allocate the q. */
    /*                                                              */
    if (l1que.nq != nq || npix != l1file->npix) {

        npix = l1file->npix;

        /* Free any previously allocated queue */
        free_l1q();

        /* Allocate a fresh queue */
        if (alloc_l1q(nq, l1file, &l1que) == 0) {
            printf("-E- %s: Unable to allocate L1B record queue.\n",
                    __FILE__);
            return (EXIT_FAILURE);
        }

        for(iq=0;iq<nq;iq++)
            l1que.r[iq].uncertainty=NULL;
        if(l1rec->uncertainty){
            for(iq=0;iq<nq;iq++){
                l1que.r[iq].uncertainty= (uncertainty_t*) malloc(sizeof (uncertainty_t));
                if (alloc_uncertainty(nbands, l1rec->uncertainty->nbands_ac, npix,l1que.r[iq].uncertainty) != 0) {
                    printf("-E- %s: Unable to allocate error record.\n",__FILE__);
                    exit(EXIT_FAILURE);
                }
            }
        }
    }

    /* Ensure that the queue is centered on iscan                   */
    if (loadl1q(l1file, iscan, dscan) != 0) {
        printf("-E- %s %d: Error reading %s at scan %d.\n",
                __FILE__, __LINE__, l1file->name, iscan);
        return (EXIT_FAILURE);
    }

    /* Now make a copy of the center queue record.                  */
    cpl1rec(l1rec, crec);

    /* And replace data with smoothed values, if desired.           */
    if (input->fctl.nfilt > 0) {
        filter(&input->fctl, &l1que, l1rec, dscan);
        setflagbits_l1(0, l1rec, -1);
        setflagbits_l1(1, l1rec, -1);
        for (ip = 0; ip < npix; ip++)
            l1_mask_set(l1rec, ip);
    }

    //lt_agregat_ocis(l1rec);

    if(l1rec->uncertainty){

        init_uncertainty(l1rec->uncertainty,1);

        /*   calculate the error in Lt and Lr */
        //set_error_input(l1rec);

        get_uncertainty(l1rec);

    }

    /* Reset scan-specific private data */
    seawater_set(l1rec);

    return (EXIT_SUCCESS);
}
