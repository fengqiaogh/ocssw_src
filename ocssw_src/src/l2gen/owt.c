#include "l12_proto.h"
#include "nr.h"
#include "nrutil.h"

#define DEBUG_OWT 0
#define NCLASSES 16
#define NRRS 6
#define NWTS 9

static float **owt;
static float **owtn;
static float *owtd;
static int nclass = NCLASSES;
static int nwts = NWTS;
static int lastScanRun = -1;

//int fuzzy_func (float *rrs, float **urrs,float ***y3inv, int nclass, int df, double *outdata);
void fuzzy_func_v3(float *rrs, float **urrs, float ***y3inv, int nclass, int nwts, int df, double *outdata);
void covariance_inversion(float *rrs_cov, int nclasses, int df, float ***y3inv);

/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
void run_owt(l2str *l2rec) {
    static int firstCall = 1;
    static int bindx[NRRS];
    static int nrrs = NRRS;
    static float badval = BAD_FLT;

    static float rrs[NRRS];
    static float rrs_means[NRRS * NCLASSES];
    static float rrs_cov[NCLASSES * NRRS * NRRS];

    static float **urrs;
    static float ***y3inv;
    static double outclass[NCLASSES];
    static float *wave;

    float owt_sum;
    float owt_max;
    int status = 0;
    int32_t ip, ib, ic, ipb;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;

    if (firstCall) {

        int32_t ic;
        int32 sd_id, sds_id, i, j, k;
        int32 dims1[2], start1[2], edge1[2];
        int32 dims2[3], start2[3], edge2[3];
        char fname[FILENAME_MAX];
        char sdsname[H4_MAX_NC_NAME] = "";

        firstCall = 0;

        strcpy(fname, input->owtfile);
        if (strlen(fname) == 0) {
            printf("-E- %s line %d: No owtfile specified.\n", __FILE__, __LINE__);
            exit(1);
        }
        printf("\nLoading optical class data from %s\n", fname);
        sd_id = SDstart(fname, DFACC_RDWR);

        // sensor specifics
        if ((wave = (float *) calloc(l1file->nbands, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to run_owt\n");
            exit(FATAL_ERROR);
        }

        switch (l1file->sensorID) {
        case SEAWIFS:
            nrrs = 5;
            wave[0] = 412;
            wave[1] = 443;
            wave[2] = 490;
            wave[3] = 510;
            wave[4] = 555;
            break;
        case MODIST:
        case MODISA:
            nrrs = 4;
            wave[0] = 412;
            wave[1] = 443;
            wave[2] = 490;
            wave[3] = 550;
            break;
        default:
            printf("%s Line %d: No classification data available for this sensor.\n", __FILE__, __LINE__);
            break;
        }

        // input band indicies
        for (ib = 0; ib < nrrs; ib++)
            bindx[ib] = windex(wave[ib], l1file->fwave, l1file->nbands);

        // static storage for results
        urrs = matrix(1, nrrs, 1, nclass);
        y3inv = f3tensor(1, nclass, 1, nrrs, 1, nrrs);
        owt = allocate2d_float(l1rec->npix, nclass);
        if (owt == NULL) {
            printf("-E- %s line %d : error allocating memory for optical classes.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        owtn = allocate2d_float(l1rec->npix, nclass);
        if (owtn == NULL) {
            printf("-E- %s line %d : error allocating memory for optical classes.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        owtd = (float *) calloc(l1rec->npix, sizeof (float));
        if (owtd == NULL) {
            printf("-E- %s line %d : error allocating memory for optical classes.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        dims1 [0] = nrrs;
        dims1 [1] = nclass;
        start1 [0] = 0;
        start1 [1] = 0;
        edge1 [0] = dims1[0];
        edge1 [1] = dims1[1];

        strcpy(sdsname, "class_means");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDreaddata(sds_id, start1, NULL, edge1, (VOIDP) rrs_means);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, fname);
            SDend(sd_id);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        dims2 [0] = nclass;
        dims2 [1] = nrrs;
        dims2 [2] = nrrs;
        start2 [0] = 0;
        start2 [1] = 0;
        start2 [2] = 0;
        edge2 [0] = dims2[0];
        edge2 [1] = dims2[1];
        edge2 [2] = dims2[2];

        strcpy(sdsname, "class_covariance");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDreaddata(sds_id, start2, NULL, edge2, (VOIDP) rrs_cov);
        if (status != 0) {
            printf("-E- %s:  Error reading SDS %s from %s.\n", __FILE__, sdsname, fname);
            SDend(sd_id);
            exit(1);
        } else {
            status = SDendaccess(sds_id);
        }

        for (ic = 0; ic < nclass; ic++)
            for (ib = 0; ib < nrrs; ib++) {
                urrs[ib + 1][ic + 1] = rrs_means[ib * nclass + ic];
            }

        covariance_inversion(rrs_cov, nclass, nrrs, y3inv);

        if (DEBUG_OWT) {
            printf("\n");
            for (i = 0; i < nclass; i++) for (j = 0; j < nrrs; j++)
                    printf("urrs %d %d %f\n", j + 1, i + 1, urrs[j + 1][i + 1]);
            printf("\n");
            for (i = 0; i < nclass; i++) for (j = 0; j < nrrs; j++) for (k = 0; k < nrrs; k++)
                        printf("y3inv %d %d %d %f\n", i + 1, j + 1, k + 1, y3inv[i + 1][j + 1][k + 1]);
        }
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        ipb = ip * l1file->nbands;

        // initialize i/o arrays

        for (ic = 0; ic < nclass; ic++) {
            owt [ip][ic] = badval;
            owtn[ip][ic] = badval;
        }
        owtd[ip] = badval;

        status = 0;
        for (ib = 0; ib < nrrs; ib++) {
            rrs[ib] = rrs_above_to_below(l2rec->Rrs[ipb + bindx[ib]]) * l2rec->Rrs[ipb + bindx[ib]];
            if (rrs[ib] < 0) status++;
            if (DEBUG_OWT) printf("\nrrs %d %f %f", ib, rrs[ib], l2rec->Rrs[ipb + bindx[ib]]);
        }

        // run classification and store
        if (status == 0) {
            // status = fuzzy_func(rrs,urrs,y3inv,nclass,nrrs,outclass);
            fuzzy_func_v3(rrs, urrs, y3inv, nclass, nwts, nrrs, outclass);
            if (status == 0) {
                owt_sum = 0.0;
                owt_max = BAD_FLT;
                for (ic = 0; ic < nwts; ic++) {
                    if (outclass[ic] < 1e-35)
                        outclass[ic] = 0.0;
                    owt[ip][ic] = (float) outclass[ic];
                    owt_sum += owt[ip][ic];
                    if (DEBUG_OWT) printf("\n%d %d %f", ip, ic, outclass[ic]);
                    if (owt[ip][ic] > owt_max) {
                        owt_max = owt[ip][ic];
                        owtd[ip] = ic + 1;
                    }
                }
                for (ic = 0; ic < nwts; ic++) {
                    if (outclass[ic] >= 0.0)
                        owtn[ip][ic] = (float) outclass[ic] / owt_sum;
                }
            }

        }

        if (status != 0)
            l1rec->flags[ip] |= PRODFAIL;
    }

    lastScanRun = l1rec->iscan;
}


/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
typedef struct error_struc {
    float rms;
    float bias;
    float perc;
} errstr;


/* ---------------------------------------------------------------------- */

/* ---------------------------------------------------------------------- */
float chl_error(char *fname, float wts[], int nwts, int dclass) {
    static int firstCall = 1;
    static errstr *erec;
    static float badval = BAD_FLT;
    static float min_sum_wts = 1e-6;

    float perc, rms, bias;
    float sum_wts;
    int ic;

    if (firstCall) {

        FILE *fp;
        char line [80];
        int clss;

        if (strlen(fname) == 0) {
            printf("-E- %s line %d: No owtchlerrfile specified.\n", __FILE__, __LINE__);
            exit(1);
        }
        printf("\nLoading chl error table from %s\n", fname);

        if ((fp = fopen(fname, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n",
                    __FILE__, __LINE__, fname);
            return (-1);
        }

        if ((erec = (errstr *) malloc(nwts * sizeof (errstr))) == NULL) {
            printf("-E- %s line %d: error allocating space for error rec.\n", __FILE__, __LINE__);
            exit(1);
        }

        ic = 0;
        while (fgets(line, 80, fp)) {
            if (line[0] == '#' || line[0] == '\n')
                continue;
            sscanf(line, "%d %f %f %f\n", &clss, &perc, &rms, &bias);
            erec[ic].perc = perc;
            erec[ic].bias = bias;
            erec[ic].rms = rms;
            ic++;
            if (ic > nwts) break;
        }

        if (ic != nwts) {
            printf("-E- %s line %d: Number of weights (%d) does not match number of error records (%d) in %s\n",
                    __FILE__, __LINE__, nwts, ic, fname);
            exit(1);
        }

        firstCall = 0;
    }

    sum_wts = 0.0;
    if (erec[dclass - 1].perc > 0.0) { // if no stats for dominant class, bail
        perc = 0.0;
        bias = 0.0;
        rms = 0.0;
        for (ic = 0; ic < nwts; ic++) {
            if (wts[ic] >= 0.0 && erec[ic].perc > 0.0) {
                perc += (wts[ic] * erec[ic].perc);
                bias += (wts[ic] * erec[ic].perc);
                rms += (wts[ic] * erec[ic].perc);
                sum_wts += wts[ic];
            }
        }
    }

    if (sum_wts > min_sum_wts)
        perc /= sum_wts;
    else
        perc = badval;

    return (perc);
}



/* ------------------------------------------------------------------- */

/* ------------------------------------------------------------------- */
void optical_water_type(l2str *l2rec, l2prodstr *p, void *vptr) {
    int classID = p->prod_ix - 1;
    int prodID = p->cat_ix;
    float *fptr = vptr;
    int32_t ip;

    l1str *l1rec = l2rec->l1rec;

    if (classID < 0 || classID > nwts) {
        printf("%s line %d: there is no class member %d\n", __FILE__, __LINE__, classID + 1);
        exit(1);
    }

    if (l1rec->iscan != lastScanRun)
        run_owt(l2rec);

    switch (prodID) {
    case CAT_owt:
        for (ip = 0; ip < l1rec->npix; ip++)
            *fptr++ = owt[ip][classID];
        break;
    case CAT_owtn:
        for (ip = 0; ip < l1rec->npix; ip++)
            *fptr++ = owtn[ip][classID];
        break;
    case CAT_owtd:
        for (ip = 0; ip < l1rec->npix; ip++)
            *fptr++ = owtd[ip];
        break;
    case CAT_chl_owterr:
        for (ip = 0; ip < l1rec->npix; ip++)
            *fptr++ = chl_error(input->owtchlerrfile, &owt[ip][0], nwts, owtd[ip]);
        break;
    }

    return;
}


