/* ============================================================================================== */
/* module water.c - functions to read and return absorption and backscatter of pure sea water     */
/*                                                                                                */
/* B. Franz, NASA/GSFC Ocean Color Discipline Processing Group, Sep. 2004                         */
/*                                                                                                */
/* ============================================================================================== */

#include "l12_proto.h"

#define MINAWTAB 200
#define MAXAWTAB 2449
#define INTAWTAB 1
#define NAWTAB  ((MAXAWTAB-MINAWTAB)/INTAWTAB + 1)

static double awtab [NAWTAB];
static double bbwtab[NAWTAB];
static int32_t ntab = NAWTAB;
static int min_wl = MINAWTAB;
static int del_wl = INTAWTAB;


/* ---------------------------------------------------------------------------------------------- */
/* read_water_spectra() - called once to load look-up table static arrays                         */
/* ---------------------------------------------------------------------------------------------- */
void read_water_spectra(void) {
    // TODO: Dynamically allocate awtab and bbwtab arrays
    static int firstCall = 1;

    FILE *fp;
    char *line;
    int32_t i, status, netcdf_input;
    float wl = 0, aw = 0, bw = 0;
    int ncid, varid;

    if (!firstCall) return;

    /* Does the file exist? */
    if (access(input->water_spectra_file, F_OK) || access(input->water_spectra_file, R_OK)) {
        printf("-E- %s: water_spectra_file '%s' does not exist or cannot open.\n",
                __FILE__, input->water_spectra_file);
        exit(EXIT_FAILURE);
    }

    /* test for NetCDF input file */
    status = nc_open(input->water_spectra_file, NC_NOWRITE, &ncid);
    netcdf_input = (status == NC_NOERR);
    if (netcdf_input) {
        /* Get the aw data. */
        if ((nc_inq_varid(ncid, "aw", &varid)) == NC_NOERR) {
            nc_get_var_double(ncid, varid, &awtab[0]);
        } else {
            printf("-E- %s: water_spectra_file '%s' does not have aw.\n",
                    __FILE__, input->water_spectra_file);
            exit(EXIT_FAILURE);
        }
        /* Get the bw data. */
        if ((nc_inq_varid(ncid, "bw", &varid)) == NC_NOERR) {
            nc_get_var_double(ncid, varid, &bbwtab[0]);
            for (i = 0; i < ntab; i++) {
                bbwtab[i] *= 0.5;
            }
        } else {
            printf("-E- %s: water_spectra_file '%s' does not have bw.\n",
                    __FILE__, input->water_spectra_file);
            exit(EXIT_FAILURE);
        }

        nc_close(ncid);
    } else {
        /* This reader is for the "old" SeaBASS formatted ASCII input
         * It still suffers from the "can't handle blank lines or missing
         * comment character" problem
         * Deprecate as soon as netCDF version is official
         */
        if ((fp = fopen(input->water_spectra_file, "r")) == NULL) {
            fprintf(stderr,
                    "-E- %s line %d: unable to open %s for reading\n", __FILE__, __LINE__, input->water_spectra_file);
            exit(1);
        }

        line = (char *) calloc(80, sizeof (char));
        i = 0;
        while (i < ntab) {
            if (fgets(line, 80, fp) == NULL) {
                fprintf(stderr, "-E- %s line %d: error reading %s at line %d\n", __FILE__, __LINE__, input->water_spectra_file, i);
                exit(1);
            }
            if (line[0] != '/' && line[0] != '!') {
                sscanf(line, "%f %f %f", &wl, &aw, &bw);
                awtab [i] = aw;
                bbwtab[i] = bw / 2.0;
                i++;
            }
        }
        free(line);
    }
    firstCall = 0;
}


/* ---------------------------------------------------------------------------------------------- */
/* aw_spectra() - returns water absorption at wavelength, wl, averaged over bandwidth, width      */

/* ---------------------------------------------------------------------------------------------- */
float aw_spectra(int32_t wl, int32_t width) {
    static int firstCall = 1;

    int32_t itab = (wl - min_wl) / del_wl;
    int32_t imin = MAX(itab - width / 2 / del_wl, 0);
    int32_t imax = MIN(itab + width / 2 / del_wl, ntab);
    float aw = 0;
    int32_t i;

    if (firstCall) {
        read_water_spectra();
        firstCall = 0;
    }

    if (itab < 0) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside aw table range.\n",__FILE__,__LINE__,wl);
        itab = 0;
    } else if (itab > ntab - 1) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside aw table range.\n",__FILE__,__LINE__,wl);
        itab = ntab - 1;
    }

    for (i = imin; i <= imax; i++)
        aw += (float) awtab[i];

    aw /= (imax - imin + 1);

    return (aw);
}


/* ---------------------------------------------------------------------------------------------- */
/* bbw_spectra() - returns water backscatter at wavelength, wl, averaged over bandwidth, width    */

/* ---------------------------------------------------------------------------------------------- */
float bbw_spectra(int32_t wl, int32_t width) {
    static int firstCall = 1;

    int32_t itab = (wl - min_wl) / del_wl;
    int32_t imin = MAX(itab - width / 2 / del_wl, 0);
    int32_t imax = MIN(itab + width / 2 / del_wl, ntab);
    float bbw = 0;
    int32_t i;

    if (firstCall) {
        read_water_spectra();
        firstCall = 0;
    }

    if (itab < 0) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside bbw table range.\n",__FILE__,__LINE__,wl);
        itab = 0;
    } else if (itab > ntab - 1) {
        //fprintf(stderr,
        //        "-W- %s line %d: wavelength of %d outside bbw table range.\n",__FILE__,__LINE__,wl);
        itab = ntab - 1;
    }

    for (i = imin; i <= imax; i++)
        bbw += (float) bbwtab[i];

    bbw /= (imax - imin + 1);

    return (bbw);
}


/* ---------------------------------------------------------------------------------------------- */
/* returns aw and bbw at specified "sensor" wavelengths, appropriate to the derived nLw           */

/* ---------------------------------------------------------------------------------------------- */
void get_aw_bbw(l2str *l2rec, float wave[], int nwave, float *aw, float *bbw) {
    int ib, iw;

    if (l1_input->outband_opt >= 2) {
        for (ib = 0; ib < nwave; ib++) {
            aw [ib] = l2rec->l1rec->sw_a[ib];
            bbw[ib] = l2rec->l1rec->sw_bb[ib];
        }
    } else {
        float *senaw = l2rec->l1rec->l1file->aw;
        float *senbbw = l2rec->l1rec->l1file->bbw;
        for (ib = 0; ib < nwave; ib++) {
            iw = bindex_get(wave[ib]);
            aw [ib] = senaw [iw];
            bbw[ib] = senbbw[iw];
        }
    }
}
