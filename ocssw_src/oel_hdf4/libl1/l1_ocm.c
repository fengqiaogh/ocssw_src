/* ============================================================================ */
/* module l1_ocm_hdf.c - functions to read OCM L1B for MSL12                    */
/*                                                                              */
/* Written By: B. Franz, NASA/OBPG, April 2008     .                            */
/*                                                                              */
/* ============================================================================ */

#include "l1_ocm.h"
#include <hdf4utils.h>

#include <hdf.h>
#include <mfhdf.h>

/* ----------------------------------------------------------------------------------- */
/* openl1_ocm() - opens a OCM L1B file for reading.                                */

/* ----------------------------------------------------------------------------------- */
int openl1_ocm(filehandle *file) {
    int32_t npix;
    int32_t nscan;
    int32_t sd_id;

    /* Open the HDF input file */
    sd_id = SDstart(file->name, DFACC_RDONLY);
    if (sd_id == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->name, DFACC_RDONLY);
        exit(1);
    }

    /* Get pixel and scan dimensions */
    if (getHDFattr(sd_id, "Number of Lines", "", (VOIDP) & nscan) != 0) {
        fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                __FILE__, __LINE__);
        exit(1);
    }
    if (getHDFattr(sd_id, "Number of Pixels", "", (VOIDP) & npix) != 0) {
        fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                __FILE__, __LINE__);
        exit(1);
    }

    file->npix = npix;
    file->ndets = 1;
    file->nscan = nscan;
    file->sd_id = sd_id;
    strcpy(file->spatialResolution, "350 m");

    return (LIFE_IS_GOOD);
}


/* ----------------------------------------------------------------------------------- */
/* readl1_ocm() - reads 1 line (scan) from an OCM L1B file, loads l1rec.           */

/* ----------------------------------------------------------------------------------- */
int readl1_ocm(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;
    static int have_ctl = 0;
    static int32_t nctlpix = 0;
    static uint16_t *data = NULL;
    static int32_t *pixnum = NULL;
    static float *ctlfpix = NULL;
    static int32_t *ctlpix = NULL;
    static float *ctllon = NULL;
    static float *ctllat = NULL;
    static float *ctlsolz = NULL;
    static float *ctlsola = NULL;
    static float *ctlsenz = NULL;
    static float *ctlsena = NULL;
    static float *dctllon = NULL;
    static float *dctllat = NULL;
    static float *dctlsolz = NULL;
    static float *dctlsola = NULL;
    static float *dctlsenz = NULL;
    static float *dctlsena = NULL;
    static int16_t factor[8];
    static char sdsname[32];
    static int16_t syear = 0;
    static int16_t sday = 0;
    static int32_t smsec = 0;
    static int32_t dmsec = 0;

    int32_t sd_id = file->sd_id;
    int32_t npix = file->npix;

    int32_t ib, ip, ipb;

    if (firstCall) {

        firstCall = 0;

        /* load control-point arrays, if required */

        if (getHDFattr(sd_id, "Number of Pixel Control Points", "", (VOIDP) & nctlpix) != 0) {
            fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if (nctlpix != npix) {
            printf("Geolocation will be derived from control points.\n");
            if (
                    (ctlfpix = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (ctlpix = (int32_t *) calloc(nctlpix, sizeof (int32_t))) == NULL ||
                    (ctllon = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (ctllat = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (ctlsolz = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (ctlsola = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (ctlsenz = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (ctlsena = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (dctllon = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (dctllat = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (dctlsolz = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (dctlsola = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (dctlsenz = (float *) calloc(nctlpix, sizeof (float))) == NULL ||
                    (dctlsena = (float *) calloc(nctlpix, sizeof (float))) == NULL) {

                printf("-E- %s line %d: Error allocating control-point space.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            READ_SDS_ID(sd_id, "ctlpixl", ctlpix, 0, 0, 0, 0, nctlpix, 1, 1, 1);
            for (ip = 0; ip < nctlpix; ip++) ctlfpix[ip] = (float) ctlpix[ip];
            have_ctl = 1;
        }

        /* get pixel numbers */
        pixnum = (int32_t *) calloc(npix, sizeof (int32_t));
        READ_SDS_ID(sd_id, "pixnum", pixnum, 0, 0, 0, 0, npix, 1, 1, 1);


        /* allocate image data */
        if ((data = (uint16_t *) calloc(npix * l1rec->l1file->nbands, sizeof (uint16_t))) == NULL) {
            printf("-E- %s line %d: Error allocating data space.\n",
                    __FILE__, __LINE__);
            return (1);
        }

        /* radiance conversion factors */
        for (ib = 0; ib < l1rec->l1file->nbands; ib++) {
            sprintf(sdsname, "band%1d", ib + 1);
            if (getHDFattr(sd_id, "factor", sdsname, (VOIDP) & factor[ib]) != 0) {
                fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                        __FILE__, __LINE__);
                exit(1);
            }
        }

        /* get time info */
        if (getHDFattr(sd_id, "Start Year", "", (VOIDP) & syear) != 0) {
            fprintf(stderr, "-E- %s line %d: error getting time info.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if (getHDFattr(sd_id, "Start Day", "", (VOIDP) & sday) != 0) {
            fprintf(stderr, "-E- %s line %d: error getting time info.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if (getHDFattr(sd_id, "Start Msec", "", (VOIDP) & smsec) != 0) {
            fprintf(stderr, "-E- %s line %d: error getting time info.\n",
                    __FILE__, __LINE__);
            exit(1);
        }
        if (getHDFattr(sd_id, "Delta Msec", "", (VOIDP) & dmsec) != 0) {
            fprintf(stderr, "-E- %s line %d: error getting time info.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

    }

    /* Get position and path geometry */
    if (have_ctl) {

        READ_SDS_ID(sd_id, "lon", ctllon, scan, 0, 0, 0, 1, nctlpix, 1, 1);
        READ_SDS_ID(sd_id, "lat", ctllat, scan, 0, 0, 0, 1, nctlpix, 1, 1);
        READ_SDS_ID(sd_id, "solz", ctlsolz, scan, 0, 0, 0, 1, nctlpix, 1, 1);
        READ_SDS_ID(sd_id, "sola", ctlsola, scan, 0, 0, 0, 1, nctlpix, 1, 1);
        READ_SDS_ID(sd_id, "senz", ctlsenz, scan, 0, 0, 0, 1, nctlpix, 1, 1);
        READ_SDS_ID(sd_id, "sena", ctlsena, scan, 0, 0, 0, 1, nctlpix, 1, 1);

        spline(ctlfpix, ctllon, nctlpix, 1e30, 1e30, dctllon);
        spline(ctlfpix, ctllat, nctlpix, 1e30, 1e30, dctllat);
        spline(ctlfpix, ctlsolz, nctlpix, 1e30, 1e30, dctlsolz);
        spline(ctlfpix, ctlsola, nctlpix, 1e30, 1e30, dctlsola);
        spline(ctlfpix, ctlsenz, nctlpix, 1e30, 1e30, dctlsenz);
        spline(ctlfpix, ctlsena, nctlpix, 1e30, 1e30, dctlsena);

        for (ip = 0; ip < npix; ip++) {
            l1rec->pixnum[ip] = pixnum[ip];
            splint(ctlfpix, ctllon, dctllon, nctlpix, (float) pixnum[ip], &l1rec->lon[ip]);
            splint(ctlfpix, ctllat, dctllat, nctlpix, (float) pixnum[ip], &l1rec->lat[ip]);
            splint(ctlfpix, ctlsolz, dctlsolz, nctlpix, (float) pixnum[ip], &l1rec->solz[ip]);
            splint(ctlfpix, ctlsola, dctlsola, nctlpix, (float) pixnum[ip], &l1rec->sola[ip]);
            splint(ctlfpix, ctlsenz, dctlsenz, nctlpix, (float) pixnum[ip], &l1rec->senz[ip]);
            splint(ctlfpix, ctlsena, dctlsena, nctlpix, (float) pixnum[ip], &l1rec->sena[ip]);
        }

    } else {
        READ_SDS_ID(sd_id, "lon", l1rec->lon, scan, 0, 0, 0, 1, npix, 1, 1);
        READ_SDS_ID(sd_id, "lat", l1rec->lat, scan, 0, 0, 0, 1, npix, 1, 1);
        READ_SDS_ID(sd_id, "solz", l1rec->solz, scan, 0, 0, 0, 1, npix, 1, 1);
        READ_SDS_ID(sd_id, "sola", l1rec->sola, scan, 0, 0, 0, 1, npix, 1, 1);
        READ_SDS_ID(sd_id, "senz", l1rec->senz, scan, 0, 0, 0, 1, npix, 1, 1);
        READ_SDS_ID(sd_id, "sena", l1rec->sena, scan, 0, 0, 0, 1, npix, 1, 1);
    }



    /* Get radiance data */
    READ_SDS_ID(sd_id, "band1", &data[0 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band2", &data[1 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band3", &data[2 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band4", &data[3 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band5", &data[4 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band6", &data[5 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band7", &data[6 * npix], scan, 0, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(sd_id, "band8", &data[7 * npix], scan, 0, 0, 0, 1, npix, 1, 1);

    for (ip = 0; ip < npix; ip++) {
        for (ib = 0; ib < l1rec->l1file->nbands; ib++) {
            ipb = ip * l1rec->l1file->nbands + ib;
            l1rec->Lt[ipb] = (float) data[ib * npix + ip] / (float) factor[ib];
            if (data[ib * npix + ip] == 0 || data[ipb] == 65535) {
                l1rec->hilt[ip] = 1;
                l1rec->Lt[ipb] = 1000.0;
            }
        }
        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;
    }

    l1rec->scantime = yds2unix(syear, sday, (double) ((smsec + scan * dmsec) / 1.e3));

    l1rec->npix = file->npix;

    return (0);
}

int closel1_ocm(filehandle *file) {
    if (SDend(file->sd_id)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
                __FILE__, __LINE__, file->sd_id, file->name);
        exit(1);
    }

    return (0);
}
