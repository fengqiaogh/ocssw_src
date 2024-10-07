#include "l12_proto.h"
#include <sys/types.h>
#include <unistd.h>
#include <nc_gridutils.h>

// module getz.c -retrieve mixed layer depth from World Ocean Atlas Climatology

static float zbad = BAD_FLT;

#define ZXANC 360
#define ZYANC 180
#define MINY   12
/* This program opens the netcdf/hdf zno3_climatology_woa1994 file located in ocssw//run/data/common. It extracts the mixed layer depth information from the hdf file. */

/* 	:Title = "WOA Nitrocline Depth Monthly Climatology"
 :Description = "World Ocean Atlas 1994, Mixed Layer Depth, Monthly"
 :Reference = "http://www.nodc.noaa.gov/OC5/WOA94/mix.html"
 */

float get_zno3(float lon, float lat, int day) {
    static int firstCall = 1;
    static int zx = ZXANC;
    static int zy = ZYANC;
    static float dx = 360.0 / ZXANC;
    static float dy = 180.0 / ZYANC;

    typedef float ref_t[ZXANC + 2];
    static ref_t* zref;
    static ref_t* z_sav;

    float zno3 = zbad;
    int i, j, jj;
    int status;
    float xx, yy;
    float t, u;

    if (firstCall) {
        float ztmp[12][ZXANC][ZYANC];

        char sdsname[H4_MAX_NC_NAME];
        char sdsname2[H4_MAX_NC_NAME];
        int ncid, ndims, nvars, ngatts, unlimdimid;
        int32 sds_id;
        nc_type rh_type; /* variable type */
        int rh_dimids[H4_MAX_VAR_DIMS]; /* dimension IDs */
        int rh_natts; /* number of attributes */
        int rh_ndims; /* number of dims */
        int32 ix, iy, ct;
        int32 lymin, lxmin, lymax, lxmax;
        float sum;

        char *OCSSW_root, zclimatology_file[300];

        firstCall = 0;

        // allocate data
        zref = (ref_t*) allocateMemory((ZYANC + 2) * sizeof (ref_t), "zref");
        z_sav = (ref_t*) allocateMemory((ZYANC + 2) * sizeof (ref_t), "z_sav");

        if ((OCSSW_root = getenv("OCDATAROOT")) == NULL) {
            printf(
                    "-E- %s:  Error looking up environmental variable OCDATAROOT\n",
                    __FILE__);
            exit(1);
        }
        strcpy(zclimatology_file, OCSSW_root);
        strcat(zclimatology_file, "/common/zno3_climatology_woa2005.nc");

        printf("\nOpening z climatology file %s\n\n", zclimatology_file);

        /* calculate the correct date and day of the  month using yd2md utility. */

        int16 mon = (int) day / 31; // month of year (no need for perfection..at least according to the sea surface salinity reference algorithm)

        /* Open the file */
        strcpy(sdsname, "zno3");

        /* try netCDF first */
        if (nc_open(zclimatology_file, NC_NOWRITE, &ncid) == NC_NOERR) {

            status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, zclimatology_file);
                handle_nc_error(status, __FILE__, __LINE__);
            }

            status = nc_inq_varid(ncid, sdsname, &sds_id);
            if (status != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, zclimatology_file);
                handle_nc_error(status, __FILE__, __LINE__);
            }
            nc_inq_varname(ncid, sds_id, sdsname2);
            status = nc_inq_var(ncid, sds_id, 0, &rh_type, &rh_ndims, rh_dimids,
                    &rh_natts);

            //printf("rh_ndims=%d rh_dimids=%d %d %d \n",rh_ndims,rh_dimids[0],rh_dimids[1],rh_dimids[2]);

            printf("\nReading z climatology file %s ndims=%d nvars=%d sds_id=%d var=%s\n\n", zclimatology_file, ndims, nvars, sds_id, sdsname2);
            /* Read the data. */
            if (nc_get_var(ncid, sds_id, &ztmp[0][0][0]) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, zclimatology_file);
                exit(1);
            }
            printf("\nClosing z climatology file %s\n\n", zclimatology_file);

            /* Close the file */
            if (nc_close(ncid) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d: error closing %s.\n",
                        __FILE__, __LINE__, zclimatology_file);
                exit(1);
            }
        } else {
            fprintf(stderr, "-E- %s line %d: error opening %s.\n",
                    __FILE__, __LINE__, zclimatology_file);
            exit(1);
        }
        /* rotate 180-deg and add wrapping border to simplify interpolation */
        /* new grid is -180.5,180.5 in i=0,361 and -90.5,90.5 in j=0,181    */


        jj = zy - 1;
        for (j = 0; j < zy; j++) {
            for (i = 0; i < zx; i++) {
                //ii = (i < zx / 2) ? i + zx / 2 : i - zx / 2;
                if (ztmp[mon][i][j] >= 0) {
                    zref[jj + 1][i + 1] = ztmp[mon][i][j];
                } else
                    zref[jj + 1][i + 1] = zbad;
                //   printf("ZREF[%d][%d]=%f\n",jj+1,i,zref[jj+1][i]);
            }
            zref[jj + 1][0] = zref[jj + 1][zx];
            zref[jj + 1][zx + 1] = zref[jj + 1][1];
            jj--;
        }
        for (i = 0; i < zx + 2; i++) {
            zref[0][i] = zref[1][i];
            zref[zy + 1][i] = zref[zy][i];
        }

        /*
         *  Extend the grid by 1 point (use only full grid boxes later)
         */



        memcpy(z_sav, zref,
                (ZYANC + 2) * (ZXANC + 2) * sizeof (float));
        for (iy = 0; iy < zy + 2; iy++) {
            lymin = (iy == 0) ? 0 : iy - 1;
            lymax = (iy == zy + 1) ? zy + 1 : iy + 1;

            for (ix = 0; ix < zx + 2; ix++) {
                if (zref[iy][ix] < zbad + 1) {
                    lxmin = (ix == 0) ? 0 : ix - 1;
                    lxmax = (ix == zx + 1) ? zx + 1 : ix + 1;

                    sum = 0.;
                    ct = 0;
                    for (j = lymin; j < lymax + 1; j++)
                        for (i = lxmin; i < lxmax + 1; i++) {
                            if (z_sav[j][i] > zbad + 1) {
                                sum += z_sav[j][i];
                                ct++;
                            }
                        }
                    if (ct > 0)
                        zref[iy][ix] = sum / ct;
                }
                //printf("zref[%d][%d]=%f monidx=%d\n",iy,ix,zref[iy][ix],mon);

            }
        }

    } /* end 1-time grid set up */


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), ZXANC + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ZYANC + 1), 0);

    /*
            //printf("i=%d j=%d %f \n",i,j,zref[j][i]);
            zno3 = zref[j][i];

            return (zno3);
     */

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* Bilinearly interpolate, only using full grid boxes */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    if ((zref[j][i] > zbad + 1) && (zref[j][i + 1] > zbad + 1)
            && (zref[j + 1][i] > zbad + 1)
            && (zref[j + 1][i + 1] > zbad + 1)) {

        zno3 = (1 - t) * (1 - u) * zref[j][i] + t * (1 - u) * zref[j][i + 1]
                + t * u * zref[j + 1][i + 1] + (1 - t) * u * zref[j + 1][i];

    } else
        zno3 = zbad;

    return (zno3);
}

