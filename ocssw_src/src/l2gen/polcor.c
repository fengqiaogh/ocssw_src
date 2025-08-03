/* ========================================================================== */
/* module polcor.c - functions to read and apply polarization correction      */
/*                                                                            */
/* Written By: B. Franz, SAIC, November 2005.                                 */
/*     W. Robinson, SAIC, 29 May 2009  generalize for use with VIIRS          */
/* ========================================================================== */

#include "l12_proto.h"
#include "polcor_hawkeye.h"
#include "polcor_oci.h"
#include <xcal.h>

#define NMIR  2  /* always 2 mirror sides                    */
#define NDET  40 /* big enough for largest detector array    */
#define NANG  15  /* big enough for largest set of table AOI or scan angles */

void polcor(l1str *l1rec, int32_t ip) {
    static int firstCall = 1;

    typedef float m_array[NMIR][NANG][NDET];
    static m_array *m12, *m13;
    static float *detfac;
    static float maxpix;
    static int32 nang;
    static float *ang = NULL, margin_s;
    static int32_t nir_s;
    static int32_t nir_l;

    static double *dm12;
    static double *dm13;
    static int *have_xcal;

    static float *pxwt; /* weights for each pixel in the line and */
    static char *pxangix; /* index of low angle storage */

    float *polcor = l1rec->polcor;
    float *delta_polcor=NULL;
    float *dpol = l1rec->dpol;

    int32_t pixnum = l1rec->pixnum[ip];
    int32_t detnum = l1rec->detnum;
    int32_t mside = l1rec->mside;
    int32_t nbands = l1rec->l1file->nbands;

    int32_t sensorID = l1rec->l1file->sensorID;

    float m1 [NANG], deriv_m1[NANG];
    float alpha;
    int32_t ib, ipb;
    int32_t idet, iang;
    float wt, zang;
    float L_x, L_qp, L_up;

    int iagsm[] = {0, 640, 1376, 3152, 4928, 5664, 6304};
    int iagpx[] = {0, 640, 1008, 1600, 2192, 2560, 3200};
    int iseq, ix1, ix2, ipx, irng, step[] = {1, 2, 3, 3, 2, 1};
    char ix, attrname[50];
    double sind = 0.0003104;
    float rad_2_deg = OEL_RADEG;
    int get_wt(float, float *, int, float *, char *);

    if(l1rec->uncertainty) {
        delta_polcor=l1rec->uncertainty->derv_pol;
    }

    if (mside < 0) mside = 0;
    if (mside > 1) mside = 1;

    for (ib = 0; ib < nbands; ib++) {
        l1rec->polcor[ip * nbands + ib] = 1.0;
    }

    if (input->pol_opt == 0 || input->pol_opt > 99)
        return;

    if (sensorID == OCI) {
        polcor_oci(l1rec, ip);
        return;
    }

    if(sensorID == HAWKEYE) {
        polcor_hawkeye(l1rec, ip);
        return;
    }

    if (firstCall) {

        char name [H4_MAX_NC_NAME] = "";
        char sdsname[H4_MAX_NC_NAME] = "";
        char file [FILENAME_MAX] = "";
        int32 sd_id;
        int32 sds_id;
        int32 rank;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[3] = {0, 0, 0};
        int32 end [3] = {1, 1, 1};
        int32 status;

        int32 ndets;
        int32 nmside;
        int32_t im, ia, id, ixb;


        if ((detfac = (float *) calloc(l1rec->l1file->nbands, sizeof (float))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(FATAL_ERROR);
        }

        if ((have_xcal = (int *) calloc(l1rec->l1file->nbands, sizeof (int))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(FATAL_ERROR);
        }

        for (ixb = 0; ixb < l1_input->xcal_nwave; ixb++) {
            if ((l1_input->xcal_opt[ixb] & XCALPOL) != 0) {
                if ((ib = bindex_get(l1_input->xcal_wave[ixb])) < 0) {
                    printf("-E- %sline %d: xcal wavelength %f does not match sensor\n",
                            __FILE__, __LINE__, l1_input->xcal_wave[ixb]);
                    exit(1);
                };
                have_xcal[ib] = 1;
            }
        }

        if ((m12 = (m_array *) calloc(l1rec->l1file->nbands, sizeof (m_array))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(FATAL_ERROR);
        }
        if ((m13 = (m_array *) calloc(l1rec->l1file->nbands, sizeof (m_array))) == NULL) {
            printf("-E- : Error allocating memory to dray_for_i\n");
            exit(FATAL_ERROR);
        }
        /* load polfile for each sensor wavelength */
        for (ib = 0; ib < nbands; ib++) {

            sprintf(file, "%s%s%d%s", input->polfile, "_", l1rec->l1file->iwave[ib], ".hdf");

            printf("Loading polarization file %s\n", file);

            /* Open the file */
            sd_id = SDstart(file, DFACC_RDONLY);
            if (sd_id == FAIL) {
                fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                        __FILE__, __LINE__, file, DFACC_RDONLY);
                exit(1);
            }

            // number of detectors in polcor file from global attr
            status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Detectors"), &ndets);
            if (status != 0) {
                printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                        __FILE__, __LINE__, "Number of Detectors", file);
                exit(1);
            }

            status = SDreadattr(sd_id, SDfindattr(sd_id, "Number of Mirror Sides"), &nmside);
            if (status != 0) {
                printf("-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                        __FILE__, __LINE__, "Number of Mirror Sides", file);
                exit(1);
            }

            if (ib == 0) {

                // decide on which angle information needs to be read
                // TODO: rework this to either use sensor name via sensorId2InstrumentName function
                // or better yet, redo the polcor files as netCDF and be consistent so this isn't required
                if ((sensorID == VIIRSN) || (sensorID == VIIRSJ1) || (sensorID == VIIRSJ2)) {
                    strcpy(attrname, "Number of Scan angles");
                    strcpy(sdsname, "scanangle");
                } else {
                    strcpy(attrname, "Number of AOIs");
                    strcpy(sdsname, "AOI");
                }

                // MODIS - attrname = Number of AOIs
                status = SDreadattr(sd_id, SDfindattr(sd_id, attrname), &nang);
                if (status != 0) {
                    printf(
                            "-E- %s Line %d:  Error reading global attribute %s from %s.\n",
                            __FILE__, __LINE__, attrname, file);
                    exit(1);
                }

                if ((ang = (float *) malloc(nang * sizeof (float))) == NULL) {
                    printf("-E- %s Line %d:  Error allocating memory for ang.\n",
                            __FILE__, __LINE__);
                    exit(1);
                }
                // MODIS - sdsname = AOI
                sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
                status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs); // grabs AOI items
                status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) ang);
                if (status != 0) {
                    printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                            __FILE__, __LINE__, sdsname, file);
                    exit(1);
                } else {
                    status = SDendaccess(sds_id);
                }
            }


            strcpy(sdsname, "m12");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            for (im = 0; im < nmside; im++) for (ia = 0; ia < nang; ia++) for (id = 0; id < ndets; id++) {
                        start[0] = im;
                        start[1] = ia;
                        start[2] = id;
                        // MODIS - read in m12  by [side][AOI][detector] save for each band ib
                        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & m12[ib][im][ia][id]);
                        if (status != 0) {
                            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                                    __FILE__, __LINE__, sdsname, file);
                            exit(1);
                        }
                    }

            status = SDendaccess(sds_id);

            strcpy(sdsname, "m13");
            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            for (im = 0; im < nmside; im++) for (ia = 0; ia < nang; ia++) for (id = 0; id < ndets; id++) {
                        start[0] = im;
                        start[1] = ia;
                        start[2] = id;
                        status = SDreaddata(sds_id, start, NULL, end, (VOIDP) & m13[ib][im][ia][id]);
                        if (status != 0) {
                            printf("-E- %s Line %d:  Error reading SDS %s from %s.\n",
                                    __FILE__, __LINE__, sdsname, file);
                            exit(1);
                        }
                    }
            status = SDendaccess(sds_id);
            status = SDend(sd_id);

            /* we need to translate from actual detectors to replicated detectors or  */
            /* aggregated detectors, to handle the HMODIS resolution switching. Get   */
            /* the conversion factor by comparing sensor ndets to table ndets.        */
            if ((sensorID == VIIRSN) || (sensorID == VIIRSJ1)  || (sensorID == VIIRSJ2))
                detfac[ib] = 1.;
            else
                detfac[ib] = l1rec->l1file->ndets * 1.0 / ndets; /* 0.25 to 4 */
        }


        /* get the max # pixels */
        if (sensorID == MODISA || sensorID == MODIST) {
            switch (l1_input->resolution) {
            case 250:
                maxpix = 5416.0;
                break;
            case 500:
                maxpix = 2708.0;
                break;
            default:
                maxpix = 1354.0;
                break;
            }
        } else if ((sensorID == VIIRSN) || (sensorID == VIIRSJ1) || (sensorID == VIIRSJ2)) {
            margin_s = l1rec->margin_s;
            if (l1rec->scn_fmt == 0)
                maxpix = 3200.;
            else
                maxpix = 6304 + margin_s * 2;
        } else {
            printf("-E- %s line %d: Mirror geometry unknown for %s\n",
                    __FILE__, __LINE__, sensorId2SensorDir(sensorID));
            exit(1);
        }

        // precompute the weights and angle (AOI on mirror or scan angle)
        // dimension index for each pixel in the scan

        if ((pxwt = (float *) malloc(maxpix * sizeof ( float))) == NULL) {
            printf("-E- %s Line %d:  Error allocating memory for weights.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if ((pxangix = (char *) malloc(maxpix * sizeof ( char))) == NULL) {
            printf("-E- %s Line %d:  Error allocating memory for angle index.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        // depending on instrument, get the angle and find the weight and index
        // into pol table for all pixels

        if ((sensorID == VIIRSN) || (sensorID == VIIRSJ1) || (sensorID == VIIRSJ2)) {
            if (l1rec->scn_fmt == 0) /* aggregated */ {
                for (irng = 0; irng < 6; irng++) {
                    ix1 = irng;
                    ix2 = irng + 1;
                    for (ipx = iagpx[ix1]; ipx < iagpx[ix2]; ipx++) {
                        iseq = ipx - iagpx[ix1];
                        zang = iagsm[ix1] + ((2 * iseq + 1) * step[irng] - 1.) / 2.;
                        zang = (zang - 3151.5) * sind * rad_2_deg;

                        // use zang in same way as below - I hate to repeat this though
                        get_wt(zang, ang, nang, &wt, &ix);
                        pxwt[ipx] = wt;
                        pxangix[ipx] = ix;
                    }
                }
            } else /* unaggregated */ {
                for (ipx = 0; ipx < maxpix; ipx++) {
                    zang = ((float) ipx - 3151.5 - margin_s) * sind * rad_2_deg;
                    get_wt(zang, ang, nang, &wt, &ix);
                    pxwt[ipx] = wt;
                    pxangix[ipx] = ix;
                }
            }

        } else {

            for (ipx = 0; ipx < maxpix; ipx++) {
                zang = 10.5 + ipx / maxpix * (65.5 - 10.5);
                get_wt(zang, ang, nang, &wt, &ix);
                pxwt[ipx] = wt;
                pxangix[ipx] = ix;
            }
        }

        nir_s = bindex_get(input->aer_wave_short);
        nir_l = bindex_get(input->aer_wave_long);

        firstCall = 0;
    }

    /* apply sensitivites to radiances to get corrections */

    for (ib = 0; ib < nbands; ib++) {

        ipb = ip * nbands + ib;
        idet = (int32_t) rint(detnum / detfac[ib]);
        alpha = l1rec->alpha[ip] / OEL_RADEG;
        L_x = l1rec->Lt[ipb] / l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb];
        if(l1rec->uncertainty) {
            delta_polcor[ipb]= 1/ l1rec->tg_sol[ipb] / l1rec->tg_sen[ipb];
        }
        if (L_x > 0.0) {
            L_qp = l1rec->L_q[ipb] * cos(2 * alpha) + l1rec->L_u[ipb] * sin(2 * alpha);
            L_up = l1rec->L_u[ipb] * cos(2 * alpha) - l1rec->L_q[ipb] * sin(2 * alpha);
            if (have_xcal[ib]) {
                dm12 = get_xcal(l1rec, XM12, l1rec->l1file->iwave[ib]);
                dm13 = get_xcal(l1rec, XM13, l1rec->l1file->iwave[ib]);
                polcor[ipb] = 1.0 / (1.0 - dm12[ip] * L_qp / L_x - dm13[ip] * L_up / L_x);
            } else {
                for (iang = pxangix[pixnum]; iang <= (pxangix[pixnum] + 1); iang++) {
                    m1[iang] = 1.0 / (1.0
                            - m12[ib][mside][iang][idet] * L_qp / L_x
                            - m13[ib][mside][iang][idet] * L_up / L_x);

                    if(l1rec->uncertainty) {
                        deriv_m1[iang]=-m1[iang]*m1[iang]*(m12[ib][mside][iang][idet] * L_qp + m13[ib][mside][iang][idet] * L_up)/(L_x*L_x);
                        deriv_m1[iang]*=delta_polcor[ipb];
                    }
                }
                polcor[ipb] = m1[(int) pxangix[pixnum]] * (1. - pxwt[pixnum]) + m1[(int) pxangix[pixnum] + 1] * pxwt[pixnum];
                if(l1rec->uncertainty) {
                    delta_polcor[ipb]=deriv_m1[(int) pxangix[pixnum]] * (1. - pxwt[pixnum]) + deriv_m1[(int) pxangix[pixnum] + 1] * pxwt[pixnum];
                }
            }
            dpol [ipb] = sqrt(pow(l1rec->L_q[ipb], 2.0) + pow(l1rec->L_u[ipb], 2.0)) / L_x;

            /* quick-fix to polcor of aerosol bands only */
            if (input->pol_opt == 6 && ib != nir_l && ib != nir_s)
                polcor[ipb] = 1.0;

        } else {
            polcor[ipb] = 1.0;
            dpol [ipb] = 0.0;
        }
    }
}

int get_wt(float zang, float *ang, int nang, float *wt, char *ix)
/*******************************************************************

   get_wt

   purpose: derive a weight and index into angle array

   Returns type: int - 0 if all is OK

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      float             zang             I      angle to interpolate to
      float *           ang              I      array of angles in 
                                                interpolation table
      int               nang             I      # of angle tie points
      float *           wt               O      derived weight
      char *            ix               O      index of first angle to use in 
                                                weighting

 *******************************************************************/ {
    int iang, ix2;

    if (zang <= ang[0])
        *ix = 0;
    else if (zang >= ang[nang - 1])
        *ix = nang - 2;
    else {
        for (iang = 1; iang < nang; iang++) {
            if (zang < ang[iang]) {
                *ix = iang - 1;
                break;
            }
        }
    }
    ix2 = *ix + 1;
    *wt = (zang - ang[(int) *ix]) / (ang[ix2] - ang[(int) *ix]);
    return 0;
}
