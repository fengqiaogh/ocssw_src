/*
 *  get_pml.c
 *
 *	MSl12 wrapper for PML IOP Algorithm
 *
 *  Plymouth Marine Laboratory, UK
 *
 * Reference: Smyth et al. (2006) "Semianalytical model for the derivation of 
 *            ocean color inherent optical properties: description, implementation,
 *            and performance assessment" Applied Optics, 45(31).
 */

#include <stdlib.h>
#include <math.h>
#include "pml.h"
#include "pml_iop.h"
#include "pml_iop_config.h"
#include "pml_iop_tables.h"
#include "pml_iop_calculate.h"

#include "l12_proto.h"
#include "l2prod.h"

static int PMLRecNum = -1;
static float badval = BAD_FLT;

static float *atot;
static float *bb;
static float *bbp;
static float *adg;
static float *aph;

static float *aw;
static float *bbw;

static int iw410;
static int iw440;
static int iw490;
static int iw510;
static int iw531;
static int iw555;
static int iw670;

int pml_ran(int recnum) {
    if (recnum == PMLRecNum)
        return 1;
    else
        return 0;
}

/*  Allocate private arrays for a single scan line */
void alloc_pml(int32_t npix, int32_t nbands) {
    atot = (float*) calloc(npix*nbands, sizeof (float));
    aph = (float*) calloc(npix*nbands, sizeof (float));
    adg = (float*) calloc(npix*nbands, sizeof (float));
    bb = (float*) calloc(npix*nbands, sizeof (float));
    bbp = (float*) calloc(npix*nbands, sizeof (float));
}

static float* alloc_bandsf(int32_t nbands, float *nbarray) {
    if ((nbarray = (float *) calloc(nbands, sizeof (float))) == NULL) {
        printf("-E- : Error allocating float memory in get_pml\n");
        exit(FATAL_ERROR);
    }
    return nbarray;
}

static double* alloc_bandsd(int32_t nbands, double *nbarray) {
    if ((nbarray = (double *) calloc(nbands, sizeof (double))) == NULL) {
        printf("-E- : Error allocating double memory in get_pml\n");
        exit(FATAL_ERROR);
    }
    return nbarray;
}

void run_pml(l2str *l2rec) {
    char *tmp_str;
    char configfname[FILENAME_MAX];

    static int firstCall = 1;
    int result = 0; /* flag for successful iterations */
    int CASEII; /* flag for caseII waters (set by the bright pixel code) */

    static float *Rrs; /* above surface Rrs */
    static double *rho_w; /* above surface water reflectance */
    static float *buf; /* temporary storage for QAA algorithm */
    static double *a_pml;
    static double *bb_pml;
    static double *bbp_pml;
    static double *adg_pml;
    static double *aph_pml;
    double solz;
    float senz, phi;

    int32_t ip, ib, ipb;

    l1str *l1rec = l2rec->l1rec;
    filehandle *l1file = l1rec->l1file;
    int32_t nbands = l1file->nbands;


    buf = alloc_bandsf(4 * nbands, buf);
    Rrs = alloc_bandsf(nbands, Rrs);
    rho_w = alloc_bandsd(nbands, rho_w);
    a_pml = alloc_bandsd(nbands, a_pml);
    bb_pml = alloc_bandsd(nbands, bb_pml);
    bbp_pml = alloc_bandsd(nbands, bbp_pml);
    adg_pml = alloc_bandsd(nbands, adg_pml);
    aph_pml = alloc_bandsd(nbands, aph_pml);

    if (firstCall) {

        firstCall = 0;
        if ((aw = (float *) calloc(nbands + 1, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to aw in get_pml\n");
            exit(FATAL_ERROR);
        }
        if ((bbw = (float *) calloc(nbands + 1, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to bbw in get_pmln");
            exit(FATAL_ERROR);
        }

        /* Load the various LUTs into memory */
        if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
            printf("OCDATAROOT environment variable is not defined.\n");
            exit(1);
        }
        strcpy(configfname, tmp_str);
        strcat(configfname, "/common/pml.cfg");
        printf("Loading PML IOP model config file: %s\n", configfname);
        load_work_tab(configfname, l1file->sensorID);

        /* Load the various parameters into memory */
        load_config(configfname);

        alloc_pml(l1rec->npix, nbands);

        iw410 = bindex_get(410);
        iw440 = bindex_get(440);
        iw490 = bindex_get(490);
        iw510 = bindex_get(510);
        iw531 = bindex_get(531);
        iw555 = bindex_get(551);
        iw670 = bindex_get(670);

        if (iw531 > 0 && iw510 < 0) {
            printf("PML IOP model using MODIS eps_a slopes\n");
            iw510 = iw531;
        }

        if (iw410 < 0 || iw440 < 0 || iw490 < 0 || iw510 < 0 || iw555 < 0) {
            printf("-E- %s line %d: PML model requires bands near 410, 440, 490, 510 and 555nm\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        if (l1_input->outband_opt >= 2) {
            for (ib = 0; ib < nbands; ib++) {
                aw [ib] = aw_spectra(l1file->iwave[ib], BANDW);
                bbw[ib] = bbw_spectra(l1file->iwave[ib], BANDW);
            }
        } else {
            for (ib = 0; ib < nbands; ib++) {
                aw [ib] = l1file->aw[ib];
                bbw[ib] = l1file->bbw[ib];
            }
        }

        if (!pml_is_initialized()) {
            pml_init(1, iw410, iw440, iw490, iw510, iw555, iw670, aw, bbw);
        }
    }

    for (ip = 0; ip < l1rec->npix; ip++) {

        /* Get the angular information */
        solz = (M_PI / 180.)*(l1rec->solz[ip]);
        senz = (M_PI / 180.)*(l1rec->senz[ip]);
        phi = (M_PI / 180.)*(l1rec->delphi[ip]);

        /* clear static globals */
        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            bb [ipb] = -1.0;
            bbp [ipb] = -1.0;
            atot[ipb] = -1.0;
            adg [ipb] = -1.0;
            aph [ipb] = -1.0;
            /* difference for these parameters as done on line-by-line basis */
            a_pml[ib] = -1.0;
            bb_pml[ib] = -1.0;
            adg_pml[ib] = -1.0;
            aph_pml[ib] = -1.0;
        }

        if (!l1rec->mask[ip]) {

            ipb = ip*nbands;

            /* compute Rrs at 412, 443, 488/490, 531/510, 551/555, 667/670 */
            for (ib = 0; ib < nbands; ib++) {
                Rrs[ib] = l2rec->Rrs[ipb + ib];
                /* convert Rrs into rho_w for running of model */
                rho_w[ib] = M_PI * Rrs[ib];
            }

            /* run the model in here */
            CASEII = 0;
            /* if (l2rec->bpsed[ip] > 0.)
               CASEII = 1; */

            result = iop_model(rho_w, solz, senz, phi, a_pml, bb_pml, adg_pml, aph_pml, iw531, CASEII);

            /* store results for this pixel in static globals */
            for (ib = 0; ib < nbands; ib++) {

                ipb = ip * nbands + ib;

                if (isfinite(bb_pml[ib]) && result != 1) {
                    bbp[ipb] = bb_pml[ib];
                    bb [ipb] = bb_pml[ib] + bbw[ib];
                } else {
                    bbp[ipb] = badval;
                    bb [ipb] = badval;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                if (isfinite(a_pml[ib]) && result != 1)
                    atot[ipb] = a_pml[ib] + aw[ib];
                else {
                    atot[ipb] = badval;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                if (isfinite(adg_pml[ib]) && result != 2)
                    adg[ipb] = adg_pml[ib];
                else {
                    adg[ipb] = badval;
                    l1rec->flags[ip] |= PRODFAIL;
                }

                if (isfinite(aph_pml[ib]) && result != 2)
                    aph[ipb] = aph_pml[ib];
                else {
                    aph[ipb] = badval;
                    l1rec->flags[ip] |= PRODFAIL;
                }
            }
        }
    }

    PMLRecNum = l1rec->iscan;

    free(buf);
    free(Rrs);
    free(rho_w);
    free(a_pml);
    free(bb_pml);
    free(bbp_pml);
    free(adg_pml);
    free(aph_pml);

    return;
}

/* Interface to l2_hdf_generic() to return PML products */

void get_pml(l2str *l2rec, l2prodstr *p, float prod[]) {
    int band = p->prod_ix;
    int prodID = p->cat_ix;
    int ip, ipb;

    if (!pml_ran(l2rec->l1rec->iscan))
        run_pml(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++) {

        ipb = ip * l2rec->l1rec->l1file->nbands + band;

        switch (prodID) {

        case CAT_a_pml:
            prod[ip] = (float) atot[ipb];
            break;

        case CAT_adg_pml:
            prod[ip] = (float) adg[ipb];
            break;

        case CAT_aph_pml:
            prod[ip] = (float) aph[ipb];
            break;

        case CAT_bb_pml:
            prod[ip] = (float) bb[ipb];
            break;

        case CAT_bbp_pml:
            prod[ip] = (float) bbp[ipb];
            break;

        default:
            printf("-E- %s line %d : erroneous product ID %d passed to get_pml().\n",
                    __FILE__, __LINE__, prodID);
            exit(1);
        }
    }

    return;
}

/* Interface to convl12() to return PML iops */
void iops_pml(l2str *l2rec) {
    int32_t ib, ip, ipb;
    int32_t nbands = l2rec->l1rec->l1file->nbands;

    if (!pml_ran(l2rec->l1rec->iscan))
        run_pml(l2rec);

    for (ip = 0; ip < l2rec->l1rec->npix; ip++)
        for (ib = 0; ib < nbands; ib++) {
            ipb = ip * nbands + ib;
            l2rec->a [ipb] = (float) atot[ipb];
            l2rec->bb[ipb] = (float) bb [ipb];
        }

    return;
}

