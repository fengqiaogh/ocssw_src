
/* ============================================================================ */
/* module l1_ocmdb_hdf.c - functions to read OCM L1B (NRL HDF) for MSL12        */
/* Written By: Paul Martinolich, Naval Research Laboratory                      */
/* ============================================================================ */

#include "l1_ocmdb.h"
#include "l1.h"
#include <hdf4utils.h>

#include <hdf.h>
#include <mfhdf.h>

static int32_t fileID;
static int32_t spix = 0;
static int32_t year, day, msec;

/* ----------------------------------------------------------------------------- */
/* openl1_ocmdb() - opens an OCM  L1B file for reading.                      */

/* ----------------------------------------------------------------------------- */

int openl1_ocmdb(filehandle *file) {
    int32_t npix;
    int32_t nscan;
    int32_t sds_id;
    int32_t rank;
    int32_t dims[3];
    int32_t type;
    int32_t numattr;
    int itmp;
    double itmp2;
    char date[16];

    /* Open the HDF input file */
    fileID = SDstart(file->name, DFACC_RDONLY);
    if (fileID == FAIL) {
        fprintf(stderr, "-E- %s line %d: SDstart(%s, %d) failed.\n",
                __FILE__, __LINE__, file->name, DFACC_RDONLY);
        return (HDF_FUNCTION_ERROR);
    }

    /* Get pixel and scan dimensions */
    sds_id = SDselect(fileID, SDnametoindex(fileID, ("ocm_ch1")));
    if (SDgetinfo(sds_id, NULL, &rank, dims, &type, &numattr) == -1) {
        fprintf(stderr, "-E- %s line %d: error getting dimension info.\n",
                __FILE__, __LINE__);
        return (HDF_FUNCTION_ERROR);
    }
    npix = dims[1];
    nscan = dims[0];

    getHDFattr(fileID, "pass_date", "", (VOIDP) & itmp);
    getHDFattr(fileID, "start_time", "", (VOIDP) & itmp2);
    sprintf(date, "%d%g", itmp, itmp2);
    printf("%s\n", date);
    date2ydmsec(date, &year, &day, &msec);

    /*
    rgb_map = get_bindx_rgb( file->input );
    if ( rgb_map )
        printf("OCM rgb bands: red %d nm (#%d), green %d nm (#%d), blue %d nm (#%d)\n",
                           file->input->rgb_wave[0], rgb_map[0],
                           file->input->rgb_wave[1], rgb_map[1],
                           file->input->rgb_wave[2], rgb_map[2]);
     */

    file->npix = npix;
    file->nscan = nscan;
    file->nbands = 8;
    file->sd_id = fileID;
    strcpy(file->spatialResolution, "350 m");

    return (LIFE_IS_GOOD);
}


/* ----------------------------------------------------------------------------- */
/* readl1_ocmdb() - reads 1 line (scan) from an OCM L1B file, loads l1rec.   */

/* ----------------------------------------------------------------------------- */
int readl1_ocmdb(filehandle *file, int32_t scan, l1str *l1rec) {
    static int firstCall = 1;

    static uint16 **dataarr;
    static float **scanarr;

    static double *m;
    static double *b;
    static uint16 **valid_range;
    static uint16 *FillValue;

    static char *names[] = {"ocm_ch1", "ocm_ch2", "ocm_ch3", "ocm_ch4",
        "ocm_ch5", "ocm_ch6", "ocm_ch7", "ocm_ch8"};

    int32_t npix = (int32_t) file->npix;
    int32_t ip, ib, ipb;

    if (firstCall) {

        firstCall = 0;

        printf("OCM has %d bands\n", l1rec->l1file->nbands);

        dataarr = (uint16 **) calloc(1, l1rec->l1file->nbands * sizeof (uint16*));
        scanarr = (float **) calloc(1, l1rec->l1file->nbands * sizeof (float*));

        m = (double *) calloc(1, l1rec->l1file->nbands * sizeof (double));
        b = (double *) calloc(1, l1rec->l1file->nbands * sizeof (double));

        valid_range = (uint16 **) calloc(1, l1rec->l1file->nbands * sizeof (uint16*));
        FillValue = (uint16 *) calloc(1, l1rec->l1file->nbands * sizeof (uint16));

        for (ib = 0; ib < l1rec->l1file->nbands; ib++) {
            if ((dataarr[ib] = (uint16 *) calloc(npix, sizeof (uint16))) == NULL) {
                printf("-E- %s line %d: Error allocating data space.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            if ((scanarr[ib] = (float *) calloc(npix, sizeof (float))) == NULL) {
                printf("-E- %s line %d: Error allocating scan space.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            if ((valid_range[ib] = (uint16 *) calloc(2, sizeof (uint16))) == NULL) {
                printf("-E- %s line %d: Error allocating scan space.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            if (getHDFattr(fileID, "scale_factor", names[ib], (VOIDP) & m[ib]) != 0) {
                printf("-E- %s line %d: Error reading reflectance scale attribute.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            if (getHDFattr(fileID, "add_offset", names[ib], (VOIDP) & b[ib]) != 0) {
                printf("-E- %s line %d: Error reading reflectance scale attribute.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            if (getHDFattr(fileID, "_FillValue", names[ib], (VOIDP) & FillValue[ib]) != 0) {
                printf("-E- %s line %d: Error reading reflectance scale attribute.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            if (getHDFattr(fileID, "valid_range", names[ib], (VOIDP) valid_range[ib]) != 0) {
                printf("-E- %s line %d: Error reading reflectance scale attribute.\n",
                        __FILE__, __LINE__);
                return (1);
            }
            printf("OCM Channel %s\n", names[ib]);
            printf("\tscale_factor %g, add_offset %g\n", m[ib], b[ib]);
            printf("\t_FillValue %d\n", (int) FillValue[ib]);
            printf("\tvalid_range [%d,%d]\n", (int) valid_range[ib][0], (int) valid_range[ib][1]);
        }

    }
    l1rec->scantime = yds2unix(year, day, (double) (msec / 1.e3));

    /* Get position and path geometry */

    READ_SDS_ID(fileID, "latitude", l1rec->lat, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(fileID, "longitude", l1rec->lon, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(fileID, "sun_zenith", l1rec->solz, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(fileID, "rel_azimuth", l1rec->sola, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(fileID, "sat_zenith", l1rec->senz, scan, spix, 0, 0, 1, npix, 1, 1);
    READ_SDS_ID(fileID, "scatter_phase", l1rec->sena, scan, spix, 0, 0, 1, npix, 1, 1);

    for (ip = 0; ip < npix; ip++) {

        l1rec->pixnum[ip] = spix + ip;

        if (l1rec->lon[ip] < -181.0 || l1rec->lon[ip] > 181.0 ||
                l1rec->lat[ip] < -91.0 || l1rec->lat[ip] > 91.0)
            l1rec->navfail[ip] = 1;
    }

    // read in the data and apply scaling

    for (ib = 0; ib < l1rec->l1file->nbands; ib++) {
        READ_SDS_ID(fileID, names[ib], &dataarr[ib][0], scan, spix, 0, 0, 1, npix, 1, 1);
        for (ip = 0; ip < npix; ip++) {
            scanarr[ib][ip] = (float) dataarr[ib][ip] * m[ib] + b[ib];
        }
    }

    /* copy over to l1rec */

    for (ip = 0; ip < npix; ip++) {
        for (ib = 0; ib < l1rec->l1file->nbands; ib++) {

            //         ipb = ip*l1rec->nbands + ib;
            ipb = ip * l1rec->l1file->nbands + ib;
            l1rec->Lt[ipb] = scanarr[ib][ip];

            if (dataarr[ib][ip] < valid_range[ib][0] || dataarr[ib][ip] > valid_range[ib][1]) {
                l1rec->hilt[ip] = 1;
                l1rec->Lt[ipb] = 1000.0;
            }
            if (dataarr[ib][ip] == FillValue[ib]) {
                l1rec->hilt[ip] = 1;
                l1rec->Lt[ipb] = 1000.0;
            }
            if (dataarr[ib][ip] == 0) {
                l1rec->hilt[ip] = 1;
                l1rec->Lt[ipb] = 1000.0;
            }
        }

        // copy over Lt_rgb bands if band-mapping defined
        /*
        if ( rgb_map )
            for (ib = 0; ib < 3; ib++ )
                l1rec->Lt_rgb[ip*3+ib] = l1rec->Lt[ip*l1rec->nbands+rgb_map[ib]];
         */
    }

    l1rec->npix = file->npix;

    return (LIFE_IS_GOOD);
}

int closel1_ocmdb(filehandle *file) {
    if (SDend(fileID)) {
        fprintf(stderr, "-E- %s line %d: SDend(%d) failed for file, %s.\n",
                __FILE__, __LINE__, fileID, file->name);
        return (HDF_FUNCTION_ERROR);
    }

    return (LIFE_IS_GOOD);
}
