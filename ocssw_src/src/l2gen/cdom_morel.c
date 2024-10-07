#include "l12_proto.h"

#define NTAB    100
#define MAXLINE 1024

static float badval = BAD_FLT;
static int32_t LastRecNum = -1;
static float adg_s = 0.018;
static float *chl;
static float *idx;


/* ---------------------------------------------------------------------------- */
/* run_cdom_morel - computes A. Morel's CDOM index and associated prodsucts.    */
/*                                                                              */
/* Reference:                                                                   */
/*                                                                              */
/* A. Morel, G. Gentili, A simple band ratio technique to quantify the colored  */
/* disolved and detrital organic material from ocean color remotely sensed data.*/
/* Rem. Sens. Env., 2009.                                                       */
/*                                                                              */
/* Implementation:                                                              */
/*                                                                              */
/*   B. Franz, NASA/OBPG, April 2009.                                           */
/*                                                                              */

/* ---------------------------------------------------------------------------- */
void run_cdom_morel(l2str *l2rec) {
    static int firstCall = 1;
    static float idxtab [NTAB][NTAB]; // CDOM index
    static float chltab [NTAB][NTAB]; // Chlorophyll
    static float xtab [NTAB]; // R412:R443 ratio
    static float ytab [NTAB]; // R490:R555 ratio
    static int nx;
    static int ny;
    static int ib412 = -1;
    static int ib443 = -1;
    static int ib490 = -1;
    static int ib555 = -1;

    float *Q0;
    float R412;
    float R443;
    float R490;
    float R555;
    float xrat, yrat;
    float t, u, w[4], wt;
    int i, j, ip, ipb;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {

        char filename[FILENAME_MAX];
        char *delim = ";";
        FILE *fp;
        char *tmp_str, *path;
        char line [MAXLINE];
        char buff [MAXLINE];

        if ((path = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }

        // number of visible bands and specific band indices

        ib412 = bindex_get(412);
        ib443 = bindex_get(443);
        ib490 = bindex_get(490);
        ib555 = bindex_get(551);
        ib555 = bindex_get(545);
        if (ib555 < 0) ib555 = bindex_get(550);
        if (ib555 < 0) ib555 = bindex_get(555);
        if (ib555 < 0) ib555 = bindex_get(560);
        if (ib412 < 0 || ib443 < 0 || ib490 < 0 || ib555 < 0) {
            printf("cdom_morel: incompatible sensor wavelengths for this algorithm\n");
            exit(1);
        }

        // read look-up table for band ratio to chl

        strcpy(filename, path);
        strcat(filename, "/common/morel_chl_R490_R555.dat");
        fp = fopen(filename, "r");
        if (!fp) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        printf("\nLoading Morel CDOM table from file %s\n", filename);

        j = 0;
        while (fgets(line, MAXLINE, fp) != NULL) {
            strcpy(buff, line);
            tmp_str = strtok(buff, delim);
            if (strcmp(tmp_str, "R412/R443") == 0) {
                i = 0;
                while ((tmp_str = strtok(NULL, delim)) != NULL) {
                    xtab[i] = atof(tmp_str);
                    i++;
                }
                nx = i;
            } else if (strcmp(tmp_str, "R490/R555") == 0) {
                ;
            } else {
                ytab[j] = atof(tmp_str);
                i = 0;
                while ((tmp_str = strtok(NULL, delim)) != NULL) {
                    if (strcmp(tmp_str, "NA") != 0)
                        chltab[i][j] = atof(tmp_str);
                    else
                        chltab[i][j] = -999;
                    i++;
                }
                if (i != nx) {
                    printf("-E- %s line %d: error reading (%s) file", __FILE__, __LINE__, filename);
                    exit(1);
                }
                j++;
            }
        }
        ny = j;

        printf("Read %d x %d entries.\n", nx, ny);

        // read look-up table for band ratio to CDOM index

        strcpy(filename, path);
        strcat(filename, "/common/morel_cdm_index.dat");
        fp = fopen(filename, "r");
        if (!fp) {
            printf("-E- %s line %d: error opening (%s) file", __FILE__, __LINE__, filename);
            exit(1);
        }

        printf("\nLoading Morel CDOM table from file %s\n", filename);

        j = 0;
        while (fgets(line, MAXLINE, fp) != NULL) {
            strcpy(buff, line);
            tmp_str = strtok(buff, delim);
            if (strcmp(tmp_str, "R412/R443") == 0) {
                i = 0;
                while ((tmp_str = strtok(NULL, delim)) != NULL) {
                    xtab[i] = atof(tmp_str);
                    i++;
                }
                nx = i;
            } else if (strcmp(tmp_str, "R490/R555") == 0) {
                ;
            } else {
                ytab[j] = atof(tmp_str);
                i = 0;
                while ((tmp_str = strtok(NULL, delim)) != NULL) {
                    if (strcmp(tmp_str, "NA") != 0)
                        idxtab[i][j] = atof(tmp_str);
                    else
                        idxtab[i][j] = -999;
                    i++;
                }
                if (i != nx) {
                    printf("-E- %s line %d: error reading (%s) filename", __FILE__, __LINE__, filename);
                    exit(1);
                }
                j++;
            }
        }
        ny = j;

        printf("Read %d x %d entries.\n", nx, ny);

        // allocate space for static arrays

        if ((chl = calloc(l1rec->npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Morel CDOM.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if ((idx = calloc(l1rec->npix, sizeof (float))) == NULL) {
            printf("-E- %s line %d : error allocating memory for Morel CDOM.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        firstCall = 0;
    }

    if ((Q0 = calloc(l1file->nbands, sizeof (float))) == NULL) {
        printf("-E- %s line %d : error allocating memory for Q0 Morel CDOM.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    // table look-up and interpolation for CDOM index and chl

    for (ip = 0; ip < l1rec->npix; ip++) {
        ipb = ip * l1file->nbands;
        idx[ip] = badval;
        chl[ip] = badval;
        if (l2rec->chl[ip] > 0.0) {
            // get Rrs and convert to irradiance reflectance
            qint_morel(l1file->fwave, l1file->nbands, 0.0, l2rec->chl[ip], Q0);
            R412 = l2rec->Rrs[ipb + ib412] * Q0[ib412];
            R443 = l2rec->Rrs[ipb + ib443] * Q0[ib443];
            R490 = l2rec->Rrs[ipb + ib490] * Q0[ib490];
            R555 = conv_rrs_to_555(l2rec->Rrs[ipb + ib555], l1file->fwave[ib555], -99,NULL) * Q0[ib555];
            if (R443 > 0.0 && R555 > 0.0) {
                // compute ratios and locate bounding table entries
                xrat = MAX(MIN(R412 / R443, xtab[nx - 1]), xtab[0]);
                yrat = MAX(MIN(R490 / R555, ytab[ny - 1]), ytab[0]);
                for (i = 0; i < nx - 2; i++)
                    if (xtab[i] > xrat)
                        break;
                for (j = 0; j < ny - 2; j++)
                    if (ytab[j] > yrat)
                        break;
                // bilinearly interpolate, weigh missing table values to zero
                t = (xrat - xtab[i]) / (xtab[i + 1] - xtab[i]);
                u = (yrat - ytab[j]) / (ytab[j + 1] - ytab[j]);

                w[0] = (idxtab[i ][j ] > -1) * (1 - t)*(1 - u);
                w[1] = (idxtab[i ][j + 1] > -1) * t * (1 - u);
                w[2] = (idxtab[i + 1][j + 1] > -1) * t*u;
                w[3] = (idxtab[i + 1][j ] > -1) * (1 - t) * u;

                wt = w[0] + w[1] + w[2] + w[3];

                if (wt > 0.0) {
                    idx[ip] = (idxtab[i ][j ] * w[0]
                            + idxtab[i ][j + 1] * w[1]
                            + idxtab[i + 1][j + 1] * w[2]
                            + idxtab[i + 1][j ] * w[3]) / wt;
                    chl[ip] = (chltab[i ][j ] * w[0]
                            + chltab[i ][j + 1] * w[1]
                            + chltab[i + 1][j + 1] * w[2]
                            + chltab[i + 1][j ] * w[3]) / wt;
                } else {
                    l1rec->flags[ip] |= PRODFAIL; // missing table values
                }
            } else {
                l1rec->flags[ip] |= PRODFAIL; // bad input Rrs
            }
        } else {
            l1rec->flags[ip] |= PRODFAIL; // bad input chl
        }
    }

    LastRecNum = l1rec->iscan;
    free(Q0);
}


/* ------------------------------------------------------------------- */
/* test if this line has been processed                                */

/* ------------------------------------------------------------------- */
int cdom_morel_ran(int recnum) {
    if (recnum == LastRecNum)
        return 1;
    else
        return 0;
}


/* ------------------------------------------------------------------- */
/* compute absorption (adg) for given chl and CDOM index               */

/* ------------------------------------------------------------------- */
float adg_morel(float chl, float idx, float wave) {
    if (chl > 0.0 && idx > badval + 1)
        return (idx * 0.065 * pow(chl, 0.63) * exp(-adg_s * (wave - 400)));
    else
        return (badval);
}


/* ------------------------------------------------------------------- */
/* compute percent CDOM for given chl and CDOM index                   */

/* ------------------------------------------------------------------- */
float pcdom_morel(float chl, float idx) {
    float adg;
    float aph;

    if (chl > 0.0 && idx > badval + 1) {
        aph = 0.03782 * pow(chl, 0.635);
        adg = adg_morel(chl, idx, 440.0);
        return (100 * adg / (adg + aph));
    } else
        return (badval);
}


/* ------------------------------------------------------------------- */
/* compute cdom-corrected chlorophyll                                  */

/* ------------------------------------------------------------------- */
float chl_cdomcorr_morel(float chl, float idx) {
    float A[] = {-73.65, -35.92, 15.3, 14.8};
    float X, Y;

    if (chl > 0.0 && idx > 0.0) {
        idx = MIN(MAX(idx, 0.1), 10.0);
        X = log10(idx);
        Y = X * (A[0] + X * (A[1] + X * (A[2] + X * A[3])));
        return (chl * (1 + Y / 100));
    } else
        return (badval);
}


/* ------------------------------------------------------------------- */
/* interface to convl12()                                              */

/* ------------------------------------------------------------------- */
void get_cdom_morel(l2str *l2rec, l2prodstr *p, float prod[]) {
    int prodID = p->cat_ix;
    int ib = p->prod_ix;

    int32_t ip;

    if (!cdom_morel_ran(l2rec->l1rec->iscan))
        run_cdom_morel(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        switch (prodID) {

        case CAT_iCDOM_morel:
            prod[ip] = (float) idx[ip];
            break;

        case CAT_chl_morel:
            prod[ip] = (float) chl[ip];
            break;

        case CAT_adg_morel:
            prod[ip] = adg_morel(chl[ip], idx[ip], l2rec->l1rec->l1file->fwave[ib]);
            break;

        case CAT_pCDOM_morel:
            prod[ip] = pcdom_morel(chl[ip], idx[ip]);
            break;

        case CAT_chl_cdomcorr_morel:
            prod[ip] = chl_cdomcorr_morel(l2rec->chl[ip], idx[ip]);
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to CDOM morel.\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        }
    }

    return;
}


