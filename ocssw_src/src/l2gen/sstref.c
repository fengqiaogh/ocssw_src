#include <sys/types.h>
#include <unistd.h>

/* ============================================================================ */
/* module sstref.c - retrieve sst reference temperature from ancillary file     */
/*                                                                              */
/* Written By: B. Franz, NASA/SIMBIOS, August 2003.                             */
/*                                                                              */
/* ============================================================================ */

/*
 * hardcoded the minatsrcnt value to 1 so I could eliminate the input parameter.
 * Only used for get_atsr anyway, and thus only valid for 2006-2009...
 */
#include "l12_proto.h"
#include <netcdf.h>
#include "l12_parms.h"
#include "l1_aci.h"
#define PATHCLIM 1
#define OISSTBIN 2
#define OISSTV2D 3
#define AMSRE3DAY 4
#define AMSREDAY 5
#define ATSR 6
#define NTEV2 7
#define AMSRE3DN 8
#define ATSRDAY 9
#define WINDSAT3DAY 10
#define WINDSATDAY 11
#define WINDSAT3DN 12
#define CMC 13
#define L3M 14

#define CtoK         273.15    /* conversion between degree C and Kelvin*/
static float sstbad = BAD_FLT;
static int32_t format = -1;

static int MiddleOfMonth[2][12] = {
    {15, 45, 74, 105, 135, 166, 196, 227, 258, 288, 319, 349},
    {15, 45, 75, 106, 136, 167, 197, 228, 259, 289, 320, 350}};


/* ----------------------------------------------------------------------------------- */
/* get_oisstv2d() - read and interpolate Reynolds 0.25-deg daily V2 netcdf OI files             */
/*                                                                                     */
/* B. Franz, SAIC, July 2008.                                                          */
/* ----------------------------------------------------------------------------------- */

#define OI4NX 1440
#define OI4NY 720

float get_oisstv2d(char *sstfile, float lon, float lat) {
    static int firstCall = 1;
    static int nx = OI4NX;
    static int ny = OI4NY;
    static float dx = 360.0 / OI4NX;
    static float dy = 180.0 / OI4NY;

    typedef float sstref_t[OI4NX + 2];
    static sstref_t* sstref;
    //static float sstref[OI4NY+2][OI4NX+2];

    float sst = sstbad;
    int i, j, ii;
    int ntmp;
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {
        sstref = (sstref_t*) malloc((OI4NY + 2) * sizeof (sstref_t));

        typedef int16 ssttmp_t[OI4NX];
        ssttmp_t *ssttmp;
        ssttmp = (ssttmp_t*) malloc(OI4NY * sizeof (ssttmp_t));

        char name [H4_MAX_NC_NAME] = "";
        char sdsname[H4_MAX_NC_NAME] = "";
        int ncid, ndims, nvars, ngatts, unlimdimid;
        int32 sd_id;
        int32 sds_id;
        int32 rank;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[4] = {0, 0, 0, 0};
        int32 status;
        float slope;
        float offset;

        firstCall = 0;

        printf("Loading Daily V2 0.25-deg OI Reynolds SST reference from %s\n", sstfile);
        printf("\n");

        /* Open the file */
        strcpy(sdsname, "sst");

        /* Try HDF first... */

        // SDstart seg faults on the mac, so need to protect the 
        // call with Hishdf()
        if (Hishdf(sstfile))
            sd_id = SDstart(sstfile, DFACC_RDONLY);
        else
            sd_id = FAIL;

        if (sd_id != FAIL) {

            sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
            status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
            status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) & ssttmp[0][0]);
            if (status != 0) {
                fprintf(stderr, "-E- %s Line %d:  Error reading SDS %s from %s.\n",
                        __FILE__, __LINE__, sdsname, sstfile);
                exit(1);
            }
            if (getHDFattr(sd_id, "scale_factor", sdsname, (VOIDP) & slope) != 0) {
                fprintf(stderr, "-E- %s line %d: error reading scale factor.\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            if (getHDFattr(sd_id, "add_offset", sdsname, (VOIDP) & offset) != 0) {
                fprintf(stderr, "-E- %s line %d: error reading scale offset.\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            status = SDendaccess(sds_id);
            status = SDend(sd_id);
        } else {
            
            /* try netCDF */
            if (nc_open(sstfile, NC_NOWRITE, &ncid) == NC_NOERR) {

                status = nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
                status = nc_inq_varid(ncid, sdsname, &sds_id);

                /* Read the data. */
                if (nc_get_var(ncid, sds_id, &ssttmp[0][0]) != NC_NOERR) {
                    fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                            __FILE__, __LINE__, sdsname, sstfile);
                    exit(1);
                }

                if (nc_get_att_float(ncid, sds_id, "scale_factor", &slope) != NC_NOERR) {
                    fprintf(stderr, "-E- %s line %d: error reading scale factor.\n",
                            __FILE__, __LINE__);
                    exit(1);
                }

                if (nc_get_att_float(ncid, sds_id, "add_offset", &offset) != NC_NOERR) {
                    fprintf(stderr, "-E- %s line %d: error reading scale offset.\n",
                            __FILE__, __LINE__);
                    exit(1);
                }
                /* Close the file */
                if (nc_close(ncid) != NC_NOERR) {
                    fprintf(stderr, "-E- %s line %d: error closing %s.\n",
                            __FILE__, __LINE__, sstfile);
                    exit(1);
                }

            } else {
                fprintf(stderr, "-E- %s line %d: error reading %s.\n",
                        __FILE__, __LINE__, sstfile);
                exit(1);
            }
        }

        /* rotate 180-deg and add wrapping border to simplify interpolation */
        /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721    */

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                ii = (i < nx / 2) ? i + nx / 2 : i - nx / 2;
                if (ssttmp[j][i] > -999)
                    sstref[j + 1][ii + 1] = ssttmp[j][i] * slope + offset;
                else
                    sstref[j + 1][ii + 1] = sstbad;
            }
            sstref[j + 1][0] = sstref[j + 1][nx];
            sstref[j + 1][nx + 1] = sstref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sstref[0] [i] = sstref[1][i];
            sstref[ny + 1][i] = sstref[ny][i];
        }

        free(ssttmp);
    }


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), OI4NX + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), OI4NY + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (sstref[j ][i ] > sstbad + 1) {
        ftmp += sstref[j ][i ];
        ntmp++;
    }
    if (sstref[j ][i + 1] > sstbad + 1) {
        ftmp += sstref[j ][i + 1];
        ntmp++;
    }
    if (sstref[j + 1][i + 1] > sstbad + 1) {
        ftmp += sstref[j + 1][i + 1];
        ntmp++;
    }
    if (sstref[j + 1][i ] > sstbad + 1) {
        ftmp += sstref[j + 1][i ];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sstref[j ][i ] > sstbad + 1 ? sstref[j ][i ] : ftmp);
        reftmp[0][1] = (sstref[j ][i + 1] > sstbad + 1 ? sstref[j ][i + 1] : ftmp);
        reftmp[1][1] = (sstref[j + 1][i + 1] > sstbad + 1 ? sstref[j + 1][i + 1] : ftmp);
        reftmp[1][0] = (sstref[j + 1][i ] > sstbad + 1 ? sstref[j + 1][i ] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        sst = sstbad;

    return (sst);
}

/* ----------------------------------------------------------------------------------- */
/* get_cmcsst() - read and interpolate CMC SST netcdf files v2 and v3                  */
/* ----------------------------------------------------------------------------------- */
float get_cmcsst(char *sstfile, float lon, float lat, int isKelvin) {
    static int firstCall = 1;
    static size_t nx = 0;
    static size_t ny = 0;
    static float dx = 0.1;
    static float dy = 0.1;
    
    static float **sstref;
    static float *latitudes;
    static float *longitudes;

    float sst = sstbad;
    int i, j;
    int ntmp;
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {
        short **ssttmp;

        char sdsname[H4_MAX_NC_NAME] = "";
        int32_t ncid, ndims, nvars, ngatts, unlimdimid;
        int32_t latDimID, lonDimID;
        int32 sds_id;
        float slope;
        float offset;

        firstCall = 0;

        printf("Loading SST reference from %s\n", sstfile);
        printf("\n");

        /* Open the file */
        if (nc_open(sstfile, NC_NOWRITE, &ncid) == NC_NOERR) {
            nc_inq(ncid, &ndims, &nvars, &ngatts, &unlimdimid);
            
            strcpy(sdsname, "analysed_sst");
            if (nc_inq_varid(ncid, sdsname, &sds_id) != NC_NOERR){
                strcpy(sdsname, "sst");
                if (nc_inq_varid(ncid, sdsname, &sds_id) != NC_NOERR){
                    printf("Whoops! something is wrong reading the SST reference file: %s\n", sstfile);
                    exit(EXIT_FAILURE);
                }
            }

            if (nc_inq_dimid(ncid, "lat", &latDimID) != NC_NOERR) {
                printf("Whoops! something is wrong reading the SST reference file: %s\n", sstfile);
                exit(EXIT_FAILURE);
            }
            nc_inq_dimlen(ncid, latDimID, &ny);
            dy = 180.0 / ny;
            if (nc_inq_dimid(ncid, "lon", &lonDimID) != NC_NOERR) {
                printf("Whoops! something is wrong reading the SST reference file: %s\n", sstfile);
                exit(EXIT_FAILURE);
            }
            nc_inq_dimlen(ncid, lonDimID, &nx);            
            dx = 360.0 / nx;
            
            /* size the arrays */
            ssttmp = allocate2d_short(ny, nx);
            sstref = allocate2d_float(ny+2, nx+2);
            
            latitudes = (float *) calloc(ny,sizeof(float));
            longitudes = (float *) calloc(nx,sizeof(float));
            
            /* Read the data. */
            if (nc_get_var_short(ncid, sds_id, &ssttmp[0][0]) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, sstfile);
                exit(1);
            }

            if (nc_get_att_float(ncid, sds_id, "scale_factor", &slope) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d: error reading scale factor.\n",
                        __FILE__, __LINE__);
                exit(1);
            }

            if (nc_get_att_float(ncid, sds_id, "add_offset", &offset) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d: error reading scale offset.\n",
                        __FILE__, __LINE__);
                exit(1);
            }
            strcpy(sdsname, "lat");
            nc_inq_varid(ncid, sdsname, &sds_id);
            if (nc_get_var_float(ncid, sds_id, latitudes) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, sstfile);
                exit(1);
            }
            strcpy(sdsname, "lon");
            nc_inq_varid(ncid, sdsname, &sds_id);
            if (nc_get_var_float(ncid, sds_id, longitudes) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d:  Error reading %s from %s.\n",
                        __FILE__, __LINE__, sdsname, sstfile);
                exit(1);
            }
            /* Close the file */
            if (nc_close(ncid) != NC_NOERR) {
                fprintf(stderr, "-E- %s line %d: error closing %s.\n",
                        __FILE__, __LINE__, sstfile);
                exit(1);
            }
        } else {
            fprintf(stderr, "-E- %s line %d: error opening %s.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }

        for (j = ny-1; j >= 0; j--) {
            for (i = 0; i < nx; i++) {
                if (ssttmp[j][i] > -999)
                    sstref[j + 1][i + 1] = ssttmp[j][i] * slope + offset;
                else
                    sstref[j + 1][i + 1] = sstbad;
            }
            sstref[j + 1][0] = sstref[j + 1][nx];
            sstref[j + 1][nx + 1] = sstref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sstref[0] [i] = sstref[1][i];
            sstref[ny + 1][i] = sstref[ny][i];
        }

        free2d_short(ssttmp);
    }

    /* locate position within reference grid */
    float sdiff = 999.0;
    int lonidx = 0;
    for (i=0; i < nx; i++){
        float cdiff = fabs(longitudes[i] - lon);
        if (cdiff > sdiff)
                break;
        if (cdiff < sdiff){
            sdiff = cdiff;
            lonidx = i;
        }     
    }
    sdiff = 999.0;
    int latidx = 0;
    for (i=0; i < ny; i++){
        float cdiff = fabs(latitudes[i] - lat);
        if (cdiff > sdiff)
            break;
        if (cdiff < sdiff){
            sdiff = cdiff;
            latidx = i;
        }     
    }

    /* compute longitude and latitude of that grid element */
    xx = longitudes[lonidx];
    yy = latitudes[latidx];

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (sstref[latidx ][lonidx ] > sstbad + 1) {
        ftmp += sstref[latidx ][lonidx ];
        ntmp++;
    }
    if (sstref[latidx][lonidx + 1] > sstbad + 1) {
        ftmp += sstref[latidx ][lonidx + 1];
        ntmp++;
    }
    if (sstref[latidx + 1][lonidx + 1] > sstbad + 1) {
        ftmp += sstref[latidx + 1][lonidx + 1];
        ntmp++;
    }
    if (sstref[latidx + 1][lonidx ] > sstbad + 1) {
        ftmp += sstref[latidx + 1][lonidx ];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sstref[latidx ][lonidx ] > sstbad + 1 ? sstref[latidx ][lonidx ] : ftmp);
        reftmp[0][1] = (sstref[latidx ][lonidx + 1] > sstbad + 1 ? sstref[latidx ][lonidx + 1] : ftmp);
        reftmp[1][1] = (sstref[latidx + 1][lonidx + 1] > sstbad + 1 ? sstref[latidx + 1][lonidx + 1] : ftmp);
        reftmp[1][0] = (sstref[latidx + 1][lonidx ] > sstbad + 1 ? sstref[latidx + 1][lonidx ] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];
        if (isKelvin)
            sst -= CtoK;

    } else
        sst = sstbad;

    return (sst);
}


/* ----------------------------------------------------------------------------------- */
/* get_oisst() - read and interpolate flat binary Reynolds OI SST                    */
/*                                                                                     */
/* The files were written in IEEE binary (big-endian). Each file contains four FORTRAN */
/* records described as follows:                                                       */
/*                                                                                     */
/*    rec 1: date and version number        (8 4-byte integer words)                   */
/*    rec 2: gridded sst values in degC     (360*180 4-byte real words)                */
/*    rec 3: normalized error variance      (360*180 4-byte real words)                */
/*    rec 4: gridded ice concentration      (360*180 1-byte integer words)             */
/*                                                                                     */
/* B. Franz, SAIC, May 2004.                                                           */
/* ----------------------------------------------------------------------------------- */

#define OINX 360
#define OINY 180

float get_oisst(char *sstfile, float lon, float lat) {
    static int firstCall = 1;
    static int nx = OINX;
    static int ny = OINY;
    static float dx = 360.0 / OINX;
    static float dy = 180.0 / OINY;

    typedef float sstref_t[OINX + 2];
    static sstref_t* sstref;
    //static float sstref[OINY+2][OINX+2];

    float sst = BAD_FLT;
    int i, j, ii;
    float xx, yy;
    float t, u;

    if (firstCall) {
        sstref = (sstref_t*) malloc((OINY + 2) * sizeof (sstref_t));

        typedef float ssttmp_t[OINX];
        ssttmp_t *ssttmp = (ssttmp_t*) malloc(OINY * sizeof (ssttmp_t));
        //float ssttmp[OINY][OINX];

        FILE *fp = NULL;

        int32_t syear, smon, sday;
        int32_t eyear, emon, eday;
        int32_t ndays, version;

        firstCall = 0;

        if ((fp = fopen(sstfile, "r")) == NULL) {
            printf("Error opening SST reference file %s for reading.\n", sstfile);
            exit(1);
        }

        if (fseek(fp, 4, SEEK_SET) < 0) {
            printf("Error reading SST reference file %s.\n", sstfile);
            exit(1);
        }
        if (fread(&syear, sizeof (int32_t), 1, fp) != 1) {
            printf("Error reading SST reference file %s.\n", sstfile);
            exit(1);
        }
        fread(&smon, sizeof (int32_t), 1, fp);
        fread(&sday, sizeof (int32_t), 1, fp);
        fread(&eyear, sizeof (int32_t), 1, fp);
        fread(&emon, sizeof (int32_t), 1, fp);
        fread(&eday, sizeof (int32_t), 1, fp);
        fread(&ndays, sizeof (int32_t), 1, fp);
        fread(&version, sizeof (int32_t), 1, fp);
        fseek(fp, 4, SEEK_CUR);

        if (endianess() == 1) {
            swapc_bytes((char *) &syear, sizeof (int32_t), 1);
            swapc_bytes((char *) &smon, sizeof (int32_t), 1);
            swapc_bytes((char *) &sday, sizeof (int32_t), 1);
            swapc_bytes((char *) &eyear, sizeof (int32_t), 1);
            swapc_bytes((char *) &emon, sizeof (int32_t), 1);
            swapc_bytes((char *) &eday, sizeof (int32_t), 1);
            swapc_bytes((char *) &ndays, sizeof (int32_t), 1);
            swapc_bytes((char *) &version, sizeof (int32_t), 1);
        }

        printf("Loading Weekly 1-degree OI Reynolds SST reference from %s\n", sstfile);
        printf("  file start date (y/m/d): %d / %d / %d\n", syear, smon, sday);
        printf("  file end   date (y/m/d): %d / %d / %d\n", eyear, emon, eday);
        printf("  days composited: %d\n", ndays);
        printf("  file version: %d\n", version);
        printf("\n");

        if (fseek(fp, 4, SEEK_CUR) < 0) {
            printf("Error reading SST reference file %s.\n", sstfile);
            exit(1);
        }
        if (fread(&ssttmp[0][0], sizeof (float), nx * ny, fp) != nx * ny) {
            printf("Error reading SST reference file %s.\n", sstfile);
            exit(1);
        }
        fclose(fp);

        if (endianess() == 1)
            swapc_bytes((char *) &ssttmp[0][0], 4, nx * ny);

        /* rotate 180-deg and add wrapping border to simplify interpolation */
        /* new grid is -180.5,180.5 in i=0,361 and -90.5,90.5 in j=0,181    */

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                ii = (i < nx / 2) ? i + nx / 2 : i - nx / 2;
                sstref[j + 1][ii + 1] = ssttmp[j][i];
            }
            sstref[j + 1][0] = sstref[j + 1][nx];
            sstref[j + 1][nx + 1] = sstref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sstref[0] [i] = sstref[1][i];
            sstref[ny + 1][i] = sstref[ny][i];
        }
        free(ssttmp);
    }


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), OINX + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), OINY + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    sst = (1 - t)*(1 - u) * sstref[j ][i ]
            + t * (1 - u) * sstref[j ][i + 1]
            + t * u * sstref[j + 1][i + 1]
            + (1 - t) * u * sstref[j + 1][i ];

    /*
    sst = sstref[j][i];
     */

    return (sst);
}

/* ----------------------------------------------------------------------------------- */
/* get_ntev2() - read and interpolate Reynolds 0.25-deg daily NTEV2 file	       */
/*                                                                                     */
/* B. Franz, SAIC, July 2008.                                                          */
/* ----------------------------------------------------------------------------------- */

#define OI4NX 1440
#define OI4NY 720

float get_ntev2(char *sstfile, float lon, float lat, int32_t year, int32_t day) {
    static int firstCall = 1;
    static int nx = OI4NX;
    static int ny = OI4NY;
    static float dx = 360.0 / OI4NX;
    static float dy = 180.0 / OI4NY;

    typedef float sstref_t[OI4NX + 2];
    static sstref_t* sstref;
    //static float sstref[OI4NY+2][OI4NX+2];
    //static float ssttmp[OI4NY][OI4NX];


    float sst = sstbad;
    int i, j, ii;
    int ntmp;
    int32_t syear, smon, sday;
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {
        sstref = (sstref_t*) malloc((OI4NY + 2) * sizeof (sstref_t));

        typedef float ssttmp_t[OI4NX];
        ssttmp_t *ssttmp;
        ssttmp = (ssttmp_t*) malloc(OI4NY * sizeof (ssttmp_t));

        /* daily file for years 2006..2009, 1440x720 */
        FILE *fp = NULL;
        firstCall = 0;

        /* just need first 1440x720 real*4 array of data from the daily file */
        /* iyr,imo,ida,refval */

        if (year < 2006 || year > 2009) {
            printf("ntev2 data only for years 2006 to 2009, not %d\n", year);
            exit(1);
        }
        printf("Loading ntev2 daily field from %s\n", sstfile);
        printf("\n");

        /* Open the file */
        if ((fp = fopen(sstfile, "r")) == NULL) {
            printf("-E- %s line %d: error opening SST reference file %s for reading.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }

        if (fseek(fp, 4, SEEK_SET) < 0) {
            printf("Error seeking SST NTEV2 reference file %s.\n", sstfile);
            exit(1);
        }
        if (fread(&syear, sizeof (int32_t), 1, fp) != 1) {
            printf("Error reading SST NTEV2 reference file %s.\n", sstfile);
            exit(1);
        }
        fread(&smon, sizeof (int32_t), 1, fp);
        fread(&sday, sizeof (int32_t), 1, fp);

        if (fread(&ssttmp[0][0], sizeof (float), nx * ny, fp) != nx * ny) {
            printf("Error reading SST NTEV2 reference file %s.\n", sstfile);
            exit(1);
        }
        fclose(fp);

        if (endianess() == 1) {
            swapc_bytes((char *) &syear, sizeof (int32_t), 1);
            swapc_bytes((char *) &smon, sizeof (int32_t), 1);
            swapc_bytes((char *) &sday, sizeof (int32_t), 1);
            swapc_bytes((char *) &ssttmp[0][0], sizeof (float), nx * ny);
        }

        printf("Loading Daily 0.25 Deg NTEV2 Reynolds SST reference from %s\n", sstfile);
        printf("  file date (y/m/d): %d / %d / %d\n", syear, smon, sday);
        printf("\n");


        /* rotate 180-deg and add wrapping border to simplify later interpolation between points */
        /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721    */

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                ii = (i < nx / 2) ? i + nx / 2 : i - nx / 2;
                sstref[j + 1][ii + 1] = ssttmp[j][i];
            }
            sstref[j + 1][0] = sstref[j + 1][nx];
            sstref[j + 1][nx + 1] = sstref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sstref[0] [i] = sstref[1][i];
            sstref[ny + 1][i] = sstref[ny][i];
        }
        free(ssttmp);
    }


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), OI4NX + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), OI4NY + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (sstref[j ][i ] > sstbad + 1) {
        ftmp += sstref[j ][i ];
        ntmp++;
    }
    if (sstref[j ][i + 1] > sstbad + 1) {
        ftmp += sstref[j ][i + 1];
        ntmp++;
    }
    if (sstref[j + 1][i + 1] > sstbad + 1) {
        ftmp += sstref[j + 1][i + 1];
        ntmp++;
    }
    if (sstref[j + 1][i ] > sstbad + 1) {
        ftmp += sstref[j + 1][i ];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sstref[j ][i ] > sstbad + 1 ? sstref[j ][i ] : ftmp);
        reftmp[0][1] = (sstref[j ][i + 1] > sstbad + 1 ? sstref[j ][i + 1] : ftmp);
        reftmp[1][1] = (sstref[j + 1][i + 1] > sstbad + 1 ? sstref[j + 1][i + 1] : ftmp);
        reftmp[1][0] = (sstref[j + 1][i ] > sstbad + 1 ? sstref[j + 1][i ] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        sst = sstbad;

    return (sst);
}

/* ----------------------------------------------------------------------------------- */
/* get_astr() - read and interpolate Reynolds 0.25-deg monthly ATSR file	   */
/*                                                                                     */
/* B. Franz, SAIC, July 2008.                                                          */
/* ----------------------------------------------------------------------------------- */

#define OI4NX 1440
#define OI4NY 720

float get_atsr(char *sstfile, float lon, float lat, int32_t year, int32_t day, int32_t minatsrcnt) {
    static int firstCall = 1;
    static int nx = OI4NX;
    static int ny = OI4NY;
    static int month;
    static float dx = 360.0 / OI4NX;
    static float dy = 180.0 / OI4NY;

    typedef float sstref_t[OI4NX + 2];
    static sstref_t* sstref;
    //static float sstref[OI4NY+2][OI4NX+2];
    //static float ssttmpb[2][OI4NY][OI4NX];
    //static float sstcnt[2][OI4NY][OI4NX];

    float sst = sstbad;
    int i, j, ii;
    int ntmp;
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {
        sstref = (sstref_t*) malloc((OI4NY + 2) * sizeof (sstref_t));

        typedef float ssttmp_t[OI4NY][OI4NX];
        ssttmp_t* ssttmpb = (ssttmp_t*) malloc(2 * sizeof (ssttmp_t));
        ssttmp_t* sstcnt = (ssttmp_t*) malloc(2 * sizeof (ssttmp_t));

        /* years 2006..2009, months 1..12, 1440x720 */
        FILE *fp = NULL;
        int otherleap;
        int endofyear;
        int numdays;
        int32_t foff[2];
        int32 status;
        int leap, day1, day2;
        int32_t reclen;
        float w1, w2;
        firstCall = 0;

        /* 4*12*3*(4+4+4+4+(4*1440*720)+4) */
        /* years 2006..2009; months 1..12; mean,std,count */
        /* 4 bytes of fortran record length; year,month,count (3 integers) */
        /* floats 1440x720 data; 4 byte fortran record end */
        /* use seek_set to find first month wanted? */
        /* use seek_cur to jump to next months mean? */

        if (year < 2006 || year > 2009) {
            printf("ATSR data only for years 2006 to 2009, not %d\n", year);
            exit(1);
        }
        printf("Loading ATSR monthly fields from %s\n", sstfile);
        printf("\n");
        leap = (isleap(year) == TRUE ? 1 : 0);
        for (month = 11; month >= 0; month--) {
            if (day > MiddleOfMonth[leap][month]) {
                break;
            }
        }
        /* each "record" is recbeg,iyr,imo,ndct,data,recend */
        /* which is 4+4+4+4+4*1440*720+4 = 4147220 bytes */
        /* years go from 2006 to 2009; months from 1 to 12 */
        /* mean is first of 3 data records so read one,skip two,read one */
        /* Jan 2006 mean is record 1, Feb 2006 mean is record 4, ... */

        /* read in two months and interpolate */
        reclen = 4147220; /* number of bytes in one record */

        /* linearly interpolate the reference values from the two months around day */
        /* w1 is the weight to use with the ref value from the earlier month */
        /* w2 is the weight to use with the ref value from the later month */
        /* The next months mean is in the next record after the previous months counts so set the 2nd offset */
        /* to skip the fortran end of record byte count, beg of rec byte count, year, month, count before the means */
        /* foff[1] will be set to -1 when day is in the first half of the first month or the */
        /*   last half of the last month and therefore there is no previous or next month to read */
        foff[1] = 0; /* last read is now counts, so next mean is next */
        if (month == -1) {
            /* day is in first half of Jan, interpolate with previous Dec */
            if (year == 2006) {
                /* There is no "previous Dec." so use first set as is */
                w1 = 1.0;
                w2 = 0.0;
                /* foff[0] is offset to means for first month whose center is before day */
                /* in this case it's the first data record so no offset */
                /* except the beg,year,mon,count that's before the data */
                foff[0] = 0;
                /* there is no previous month so set 2nd offset to be zero so no read is done */
                foff[1] = -1;
            } else {
                /* day is in first half of Jan, interpolate with previous Dec */
                otherleap = (isleap(year - 1) == TRUE ? 1 : 0);
                endofyear = (otherleap == 1 ? 366 : 365);
                numdays = ((endofyear - MiddleOfMonth[otherleap][11]) +
                        MiddleOfMonth[leap][0]);
                w1 = (float) (MiddleOfMonth[leap][0] - day) / (float) numdays;
                w2 = 1.0 - w1;
                /* foff[0] is offset to means for first month whose center is before day */
                /* in this case it's Dec of the previous year */
                /* data is after 16 bytes of beg+year+month+count */
                foff[0] = ((year - 1 - 2006)*12 + 11) * reclen * 3;
            }
        } else if (month == 11) {
            /* day is in last half of Dec, interpolate with next Jan */
            if (year == 2009) {
                /* there is no "next Jan." so just use Dec. coeffs as is */
                w1 = 1.0;
                w2 = 0.0;
                /* foff[0] is offset to means for first month whose center is before day */
                /* in this case it's Dec of the last year, so there is no next month */
                /* plus the stuff at the beginning of each record: beg,year,month,count */
                foff[0] = ((year - 2006)*12 + 11) * reclen * 3;
                foff[1] = -1;
            } else {
                /* day is in last half of Dec, interpolate with next Jan */
                otherleap = (isleap(year + 1) == TRUE ? 1 : 0);
                endofyear = (leap == 1 ? 366 : 365);
                numdays = ((endofyear - MiddleOfMonth[leap][11]) +
                        MiddleOfMonth[otherleap][0]);
                w2 = (float) (day - MiddleOfMonth[leap][11]) / (float) numdays;
                w1 = 1.0 - w2;
                /* foff[0] is offset to means for first month whose center is before day */
                /* in this case the first month to read is this Dec */
                /* plus fortran record count,year,month,count=16 bytes before data starts */
                foff[0] = ((year - 2006)*12 + 11) * reclen * 3;
            }
        } else {
            /* interpolate with next month */
            day1 = MiddleOfMonth[leap][month];
            day2 = MiddleOfMonth[leap][month + 1];
            w2 = (float) (day - day1) / (float) (day2 - day1);
            w1 = 1.0 - w2;
            /* foff[0] is offset to means for first month whose center is before day */
            /* in this case it's a "normal" month in the "middle" of a year */
            /* plus fortran record count,year,month,count=16 bytes before data starts */
            foff[0] = ((year - 2006)*12 + month) * reclen * 3;
        }
        /* Open the file */
        if ((fp = fopen(sstfile, "r")) == NULL) {
            printf("-E- %s line %d: error opening SST reference file %s for reading.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }
        /* skip to the first month (and skip fortran begin record count,year,month,and day count) */
        if (fseek(fp, foff[0] + 16, SEEK_SET) < 0) {
            printf("Error seeking SST reference file %s.\n", sstfile);
            exit(1);
        }

        if ((status = fread(&ssttmpb[0][0][0], sizeof (float), nx * ny, fp)) != nx * ny) {
            printf("Wrong atsr data read size: want %d got %d\n", nx*ny, status);
            exit(1);
        }
        /* skip fortran end of record, std. dev., fortran beg of record to read counts,year,month,day count */
        if (fseek(fp, 4 + reclen + 16, SEEK_CUR) < 0) {
            printf("Error seeking first counts from SST reference file %s.\n", sstfile);
            exit(1);
        }
        if ((status = fread(&sstcnt[0][0][0], sizeof (float), nx * ny, fp)) != nx * ny) {
            printf("Wrong atsr counts data read size: want %d got %d\n", nx*ny, status);
            exit(1);
        }
        if (foff[1] == 0) {
            /* read next months mean */
            /* skip fortran end of record count, begin of record count,year,month, and day count */
            fseek(fp, 20, SEEK_CUR);
            if ((status = fread(&ssttmpb[1][0][0], sizeof (float), nx * ny, fp)) != nx * ny) {
                printf("Wrong 2nd atsr data read size: want %d got %d\n", nx*ny, status);
                exit(1);
            }
            /* skip fortran end of record count for means, std. dev., begin of record count,year,month,day count to read counts */
            if (fseek(fp, 4 + reclen + 16, SEEK_CUR) < 0) {
                printf("Error seek to second counts from SST reference file %s.\n", sstfile);
                exit(1);
            }
            if ((status = fread(&sstcnt[1][0][0], sizeof (float), nx * ny, fp)) != nx * ny) {
                printf("Wrong 2nd atsr counts data read size: want %d got %d\n", nx*ny, status);
                exit(1);
            }
        }

        if (endianess() == 1) {
            swapc_bytes((char *) &ssttmpb[0][0][0], sizeof (float), 2 * nx * ny);
            swapc_bytes((char *) &sstcnt[0][0][0], sizeof (float), 2 * nx * ny);
        }

        fclose(fp);
        /* interpolate between the two months in ssttmpb */
        /* rotate 180-deg and add wrapping border to simplify later interpolation between points */
        /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721    */

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                ii = (i < nx / 2) ? i + nx / 2 : i - nx / 2;
                if (ssttmpb[0][j][i] > -999.0 && ssttmpb[1][j][i] > -999.0 &&
                        sstcnt[0][j][i] >= minatsrcnt && sstcnt[1][j][i] >= minatsrcnt)
                    sstref[j + 1][ii + 1] = w1 * ssttmpb[0][j][i] + w2 * ssttmpb[1][j][i];
                else {
                    if (ssttmpb[0][j][i] > -999.0 && sstcnt[0][j][i] >= minatsrcnt)
                        sstref[j + 1][ii + 1] = ssttmpb[0][j][i];
                    else if (ssttmpb[1][j][i] > -999.0 && sstcnt[1][j][i] >= minatsrcnt)
                        sstref[j + 1][ii + 1] = ssttmpb[1][j][i];
                    else
                        sstref[j + 1][ii + 1] = sstbad;
                }
            }
            sstref[j + 1][0] = sstref[j + 1][nx];
            sstref[j + 1][nx + 1] = sstref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sstref[0] [i] = sstref[1][i];
            sstref[ny + 1][i] = sstref[ny][i];
        }
        free(ssttmpb);
        free(sstcnt);
    }


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), OI4NX + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), OI4NY + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (sstref[j ][i ] > sstbad + 1) {
        ftmp += sstref[j ][i ];
        ntmp++;
    }
    if (sstref[j ][i + 1] > sstbad + 1) {
        ftmp += sstref[j ][i + 1];
        ntmp++;
    }
    if (sstref[j + 1][i + 1] > sstbad + 1) {
        ftmp += sstref[j + 1][i + 1];
        ntmp++;
    }
    if (sstref[j + 1][i ] > sstbad + 1) {
        ftmp += sstref[j + 1][i ];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sstref[j ][i ] > sstbad + 1 ? sstref[j ][i ] : ftmp);
        reftmp[0][1] = (sstref[j ][i + 1] > sstbad + 1 ? sstref[j ][i + 1] : ftmp);
        reftmp[1][1] = (sstref[j + 1][i + 1] > sstbad + 1 ? sstref[j + 1][i + 1] : ftmp);
        reftmp[1][0] = (sstref[j + 1][i ] > sstbad + 1 ? sstref[j + 1][i ] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        sst = sstbad;

    return (sst);
}

/* ----------------------------------------------------------------------------------- */
/* get_amsre() - read and interpolate AMSR-E data from SSMI                            */
/*                                                                                     */
/* S. Walsh, RSMAS, Aug 2010.                                                          */
/* ----------------------------------------------------------------------------------- */

#define OI4NX 1440
#define OI4NY 720

float get_amsre(char *sstfile, float lon, float lat, float solz, int32_t xsatid, int32_t sensorID) {
    static int firstCall = 1;
    static int nx = OI4NX;
    static int ny = OI4NY;
    static float dx = 360.0 / OI4NX;
    static float dy = 180.0 / OI4NY;

    typedef float sstref_t[OI4NX + 2][2];
    static sstref_t* sstref;
    //static float sstref[OI4NY+2][OI4NX+2][2];

    //    char dbgfile[FILENAME_MAX] = "./sstref.dbg.file";

    float sst = sstbad;
    int jj, ii, kk, ij;
    int ad;
    int kkmax;
    int ntmp;
    int32_t foff[2];
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {
        sstref = (sstref_t*) malloc((OI4NY + 2) * sizeof (sstref_t));

        FILE *fp = NULL;
        unsigned char ssttmp[OI4NX];
        int32 status;

        firstCall = 0;

        if (format == AMSRE3DAY)
            printf("Loading 3-Day 0.25-deg AMSR-E SST reference from %s\n", sstfile);
        else if (format == AMSRE3DN)
            printf("Loading 3-Day Day or Night 0.25-deg AMSR-E SST reference from %s\n", sstfile);
        else
            printf("Loading Daily 0.25-deg AMSR-E SST reference from %s\n", sstfile);
        printf("\n");

        /* Open the file */
        if ((fp = fopen(sstfile, "r")) == NULL) {
            printf("-E- %s line %d: error opening SST reference file %s for reading.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }

        if (format == AMSRE3DAY) {
            /* 3-Day - first field is 3 day avg of sst asc and dsc */
            foff[0] = 0;
            /* only one set for both ascending and descending */
            kkmax = 1;
        } else if (format == AMSRE3DN) {
            /* 3-Day - Day or Night - first field is 3 day avg of sst asc (day for aqua/amsr) and 2nd field is dsc night */
            foff[0] = 0;
            foff[1] = (int32_t) (nx * ny);
            /* one set of data for ascending, and one set for descending */
            kkmax = 2;
        } else {
            /* Daily - second field is 3 day avg of sst asc and dsc */
            /* v5 has 6 fields for day and 6 for night so fields 2 and 8 are sst */
            /* v7 has 7 fields for day and 7 for night so fields 2 and 9 are sst */
            foff[0] = (int32_t) (nx * ny);
            foff[1] = (int32_t) (nx * ny * 8); /* was 7 for v5 */
            /* one set of data for ascending, and one set for descending */
            kkmax = 2;
        }

        for (kk = 0; kk < kkmax; kk = kk + 1) {
            /* only read sst */
            /* ascending (day) is first, then descending (night) */
            if (fseek(fp, foff[kk], SEEK_SET) < 0) {
                fprintf(stderr, "-E- %s line %d: error reading SST reference file %s.\n",
                        __FILE__, __LINE__, sstfile);
                exit(1);
            }
            for (jj = 0; jj < ny; jj = jj + 1) {
                if ((status = fread(&ssttmp[0], sizeof (char), nx, fp)) != nx) {
                    printf("Wrong AMSR-E data read size: want %d got %d\n", nx, status);
                    exit(1);
                }
                /* 0-250 are valid values */
                for (ii = 0; ii < nx; ii = ii + 1) {
                    /* rotate 180-deg and add wrapping border to simplify interpolation */
                    /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721    */

                    ij = (ii < nx / 2) ? ii + nx / 2 : ii - nx / 2;
                    if (ssttmp[ii] <= 250) {
                        sstref[jj + 1][ij + 1][kk] = ssttmp[ii] * 0.15 - 3.0;
                    } else {
                        sstref[jj + 1][ij + 1][kk] = sstbad;
                    }
                    /* if d3d file put values in both asc and dsc slots */
                    if (kkmax == 1) {
                        sstref[jj + 1][ij + 1][1] = sstref[jj + 1][ij + 1][0];
                    }
                }

                /* wrap edges */
                sstref[jj + 1][0] [kk] = sstref[jj + 1][nx][kk];
                sstref[jj + 1][nx + 1][kk] = sstref[jj + 1][1] [kk];
            }
            /* copy edges for interpolation */
            for (ii = 0; ii < nx + 2; ii = ii + 1) {
                sstref[0] [ii][kk] = sstref[1] [ii][kk];
                sstref[ny + 1][ii][kk] = sstref[ny][ii][kk];
            }
        }
        fclose(fp);

    }

    if (solz < 90.0) {
        ad = 0; /* day - ascending for amsr on aqua */
    } else {
        ad = 1; /* night - descending for amsr on aqua */
    }

    // morning satellites, daytime pixels use night amsr
    // afternoon satellites, daytime pixels use day amsr

    // /* some day may want to swap day/night for satellites that
    //  * descend in the daytime so that the swath are the same direction
    //  * amsr was on aqua so day is 1:30 pm local time, afternoon, lots of daytime warming
    //  * so morning satellites, whose "day" pixels are before daytime warning,
    //  * should use amsr night for their reference field instead of amsr (afternoon) day
    // */

    if ((format == AMSREDAY || format == AMSRE3DAY || format == AMSRE3DN) &&
            ((sensorID == MODIST ||
            (sensorID == AVHRR && (xsatid == NO06 || xsatid == NO08 ||
            xsatid == NO10 || xsatid == NO12 || xsatid == NO15 || xsatid == NO17))))) {
        /* These satellites descend in day, opposite of aqua */
        /* and their "day" pixels are 9:30 am, before daytime warming, so use amsr night field */
        ad = 1 - ad;
    }

    /* locate LL position within reference grid */
    ii = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), nx + 1), 0);
    jj = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ny + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = ii * dx - 180.0 - dx / 2;
    yy = jj * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (sstref[jj ][ii ][ad] > sstbad + 1) {
        ftmp += sstref[jj ][ii ][ad];
        ntmp++;
    }
    if (sstref[jj ][ii + 1][ad] > sstbad + 1) {
        ftmp += sstref[jj ][ii + 1][ad];
        ntmp++;
    }
    if (sstref[jj + 1][ii + 1][ad] > sstbad + 1) {
        ftmp += sstref[jj + 1][ii + 1][ad];
        ntmp++;
    }
    if (sstref[jj + 1][ii ][ad] > sstbad + 1) {
        ftmp += sstref[jj + 1][ii ][ad];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sstref[jj ][ii ][ad] > sstbad + 1 ? sstref[jj ][ii ][ad] : ftmp);
        reftmp[0][1] = (sstref[jj ][ii + 1][ad] > sstbad + 1 ? sstref[jj ][ii + 1][ad] : ftmp);
        reftmp[1][1] = (sstref[jj + 1][ii + 1][ad] > sstbad + 1 ? sstref[jj + 1][ii + 1][ad] : ftmp);
        reftmp[1][0] = (sstref[jj + 1][ii ][ad] > sstbad + 1 ? sstref[jj + 1][ii ][ad] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        sst = sstbad;

    return (sst);
}

/* ----------------------------------------------------------------------------------- */
/* get_windsat() - read and interpolate WindSat data from SSMI                         */
/*                                                                                     */
/* S. Walsh, RSMAS, Aug 2010.                                                          */
/* ----------------------------------------------------------------------------------- */

#define OI4NX 1440
#define OI4NY 720

float get_windsat(char *sstfile, float lon, float lat, float solz, int32_t xsatid, int32_t sensorID) {
    static int firstCall = 1;
    static int nx = OI4NX;
    static int ny = OI4NY;
    static float dx = 360.0 / OI4NX;
    static float dy = 180.0 / OI4NY;

    typedef float sstref_t[OI4NX + 2][2];
    static sstref_t* sstref;
    //static float sstref[OI4NY+2][OI4NX+2][2];

    //    char dbgfile[FILENAME_MAX] = "./sstref.dbg.file";

    float sst = sstbad;
    int jj, ii, kk, ij;
    int ad;
    int kkmax;
    int ntmp;
    int32_t foff[2];
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {
        sstref = (sstref_t*) malloc((OI4NY + 2) * sizeof (sstref_t));

        FILE *fp = NULL;
        unsigned char ssttmp[OI4NX];
        float fssttmp[OI4NX];
        int32 status;

        firstCall = 0;

        if (format == WINDSAT3DAY)
            printf("Loading 3-Day 0.25-deg WindSat SST reference from %s\n", sstfile);
        else if (format == WINDSAT3DN)
            printf("Loading 3-Day Day or Night 0.25-deg WindSat SST reference from %s\n", sstfile);
        else
            printf("Loading Daily 0.25-deg WindSat SST reference from %s\n", sstfile);
        printf("\n");

        /* Open the file */
        if ((fp = fopen(sstfile, "r")) == NULL) {
            printf("-E- %s line %d: error opening SST reference file %s for reading.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }

        if (format == WINDSAT3DAY) {
            /* 3-Day - first field is 3 day avg of sst dsc (day) and asc (night) */
            foff[0] = 0;
            /* only one set for both descending and ascending */
            kkmax = 1;
        } else if (format == WINDSAT3DN) {
            /* 3-Day - Day or Night - first field is 3 day avg of sst dsc and 2nd field is asc night */
            /* values are floats */
            foff[0] = 0;
            foff[1] = (int32_t) (nx * ny * sizeof (float));
            /* one set of data for descending, and one set for ascending */
            kkmax = 2;
        } else {
            /* Daily - second field is 3 day avg of sst desc and asc */
            /* 9 fields for desc, 9 for asc, sst's are bands 2 and 11 (one based) */
            /* values are bytes */
            foff[0] = (int32_t) (nx * ny);
            foff[1] = (int32_t) (nx * ny * 10);
            /* seek_set starts at 0 each time, so skip fields 1..10 */
            /* one set of data for descending, and one set for ascending */
            kkmax = 2;
        }

        for (kk = 0; kk < kkmax; kk = kk + 1) {
            /* only read sst */
            /* descending (morning) is first, then ascending (evening) */
            if (fseek(fp, foff[kk], SEEK_SET) < 0) {
                fprintf(stderr, "-E- %s line %d: error reading SST reference file %s.\n",
                        __FILE__, __LINE__, sstfile);
                exit(1);
            }
            for (jj = 0; jj < ny; jj = jj + 1) {
                if (format == WINDSAT3DN) {
                    if ((status = fread(&fssttmp[0], sizeof (float), nx, fp)) != nx) {
                        printf("Wrong WINDSAT data read size: want %d got %d\n", nx, status);
                        exit(1);
                    }
                    //		    if (jj == 0) {
                    //		    printf(" first windsat values: %f %f %f %f %f %f %f\n", fssttmp[0],fssttmp[1],fssttmp[2],fssttmp[3],fssttmp[4],fssttmp[5],fssttmp[6]);
                    //		    }
                    if (endianess() == 1) {
                        swapc_bytes((char *) &fssttmp[0], sizeof (float), nx);
                    }
                    //		    if (jj == 220) {
                    //		        printf(" fssttmp[1164]=%f fssttmp[1165]=%f\n",fssttmp[1164],fssttmp[1165]);
                    //		        printf(" fssttmp[444]=%f fssttmp[445]=%f\n",fssttmp[444],fssttmp[445]);
                    //		    }
                    //		    if (jj == 0) {
                    //		    printf(" first windsat values: %f %f %f %f %f %f %f\n", fssttmp[0],fssttmp[1],fssttmp[2],fssttmp[3],fssttmp[4],fssttmp[5],fssttmp[6]);
                    //		    }
                } else {
                    if ((status = fread(&ssttmp[0], sizeof (char), nx, fp)) != nx) {
                        printf("Wrong WINDSAT data read size: want %d got %d\n", nx, status);
                        exit(1);
                    }
                }
                /* 0-250 are valid values */
                for (ii = 0; ii < nx; ii = ii + 1) {
                    /* rotate 180-deg and add wrapping border to simplify interpolation */
                    /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721    */

                    ij = (ii < nx / 2) ? ii + nx / 2 : ii - nx / 2;
                    if (format == WINDSAT3DN) {
                        /* for our 3 day average, bad values are -9999.0 */
                        if (fssttmp[ii] > -9998.0) {
                            sstref[jj + 1][ij + 1][kk] = fssttmp[ii];
                        } else {
                            sstref[jj + 1][ij + 1][kk] = sstbad;
                        }
                    } else {
                        if (ssttmp[ii] <= 250) {
                            sstref[jj + 1][ij + 1][kk] = ssttmp[ii] * 0.15 - 3.0;
                        } else {
                            sstref[jj + 1][ij + 1][kk] = sstbad;
                        }
                    }
                    /* if d3d file put values in both asc and dsc slots */
                    if (kkmax == 1) {
                        sstref[jj + 1][ij + 1][1] = sstref[jj + 1][ij + 1][0];
                    }
                }

                /* wrap edges */
                sstref[jj + 1][0] [kk] = sstref[jj + 1][nx][kk];
                sstref[jj + 1][nx + 1][kk] = sstref[jj + 1][1] [kk];
            }
            /* copy edges for interpolation */
            for (ii = 0; ii < nx + 2; ii = ii + 1) {
                sstref[0] [ii][kk] = sstref[1] [ii][kk];
                sstref[ny + 1][ii][kk] = sstref[ny][ii][kk];
            }
        }
        fclose(fp);


    }

    /* default is for afternoon satellites whose day pixels are around 1:30 pm local time */
    /* their day pixels should be compared to windsat evening pixels, and their night to windsat morning */
    if (solz < 90.0) {
        ad = 1; /* day pixel, use windsat evening - aescending */
    } else {
        ad = 0; /* night pixel, use windsat morning - descending */
    }

    // windsat fields are morning, descending (08:45) and evening, ascending (20:45)
    // viirs is descending at 13:45 and ascending at 01:45
    // so for viirs descending day, use evening ascending windsat
    // morning fields are best for night pixels from afternoon satellites because no daytime warming has occurred yet
    // evening fields are best for afternoon daytime pixels

    // morning satellites, daytime pixels use morning windsat
    // afternoon satellites, daytime pixels use evening windsat

    // /* some day may want to swap day/night for satellites that
    //  * descend in the daytime so that the swath are the same direction
    // */
    if ((format == WINDSATDAY || format == WINDSAT3DAY || format == WINDSAT3DN) &&
            ((sensorID == MODIST ||
            (sensorID == AVHRR && (xsatid == NO06 || xsatid == NO08 ||
            xsatid == NO10 || xsatid == NO12 || xsatid == NO15 || xsatid == NO17))))) {
        /* These satellites descend in day, opposite of aqua */
        /* and their "day" pixels are 9:30 am, before daytime warming, so use amsr night field */
        ad = 1 - ad;
    }

    /* locate LL position within reference grid */
    ii = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), nx + 1), 0);
    jj = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ny + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = ii * dx - 180.0 - dx / 2;
    yy = jj * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    if (sstref[jj ][ii ][ad] > sstbad + 1) {
        ftmp += sstref[jj ][ii ][ad];
        ntmp++;
    }
    if (sstref[jj ][ii + 1][ad] > sstbad + 1) {
        ftmp += sstref[jj ][ii + 1][ad];
        ntmp++;
    }
    if (sstref[jj + 1][ii + 1][ad] > sstbad + 1) {
        ftmp += sstref[jj + 1][ii + 1][ad];
        ntmp++;
    }
    if (sstref[jj + 1][ii ][ad] > sstbad + 1) {
        ftmp += sstref[jj + 1][ii ][ad];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sstref[jj ][ii ][ad] > sstbad + 1 ? sstref[jj ][ii ][ad] : ftmp);
        reftmp[0][1] = (sstref[jj ][ii + 1][ad] > sstbad + 1 ? sstref[jj ][ii + 1][ad] : ftmp);
        reftmp[1][1] = (sstref[jj + 1][ii + 1][ad] > sstbad + 1 ? sstref[jj + 1][ii + 1][ad] : ftmp);
        reftmp[1][0] = (sstref[jj + 1][ii ][ad] > sstbad + 1 ? sstref[jj + 1][ii ][ad] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        sst = sstbad;

    return (sst);
}

/* ----------------------------------------------------------------------------------- */
/* get_atsrday() - read and interpolate ATSR 0.1-deg daily netcdf files                */
/*                                                                                     */
/* S. Walsh, RSMAS, Feb 2011.                                                          */
/* ----------------------------------------------------------------------------------- */

#define OI1NX 3600
#define OI1NY 1800

float get_atsrday(char *sstfile, float lon, float lat) {
    static int firstCall = 1;
    static int nx = OI1NX;
    static int ny = OI1NY;
    static int32_t slen = 0;
    static float *sref = NULL;
    static float dx = 360.0 / OI1NX;
    static float dy = 180.0 / OI1NY;

    float sst = sstbad;
    int i, j, jl, nl;
    int ntmp;
    float xx, yy;
    float t, u;
    float reftmp[2][2];
    float ftmp;

    if (firstCall) {

        char name [H4_MAX_NC_NAME] = "";
        char sdsname[H4_MAX_NC_NAME] = "";
        int ll, jj, jjl;
        int32 sd_id;
        int32 sds_id;
        int32 rank;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[4] = {0, 0, 0, 0};
        int32 status;
        int32_t stlen = 0;
        float *stmp = NULL;

        firstCall = 0;

        printf("Loading Daily 0.1-deg ATSR SST reference from %s\n", sstfile);
        printf("\n");

        /* allocate arrays */
        slen = (OI1NY + 2) * (OI1NX + 2);
        if ((sref = (float *) malloc(slen * sizeof (float))) == NULL) {
            fprintf(stderr, "-E- %s line %d: Unable to allocate sstref array.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        stlen = OI1NY * OI1NX;
        if ((stmp = (float *) malloc(stlen * sizeof (float))) == NULL) {
            fprintf(stderr, "-E- %s line %d: Unable to allocate ssttmp array.\n",
                    __FILE__, __LINE__);
            exit(1);
        }

        /* Open the file */
        sd_id = SDstart(sstfile, DFACC_RDONLY);
        if (sd_id == FAIL) {
            fprintf(stderr, "-E- %s line %d: error reading %s.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }

        strcpy(sdsname, "sst_skin");
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) stmp);
        if (status != 0) {
            fprintf(stderr, "-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, sstfile);
            exit(1);
        }
        status = SDendaccess(sds_id);
        status = SDend(sd_id);

        /* ATSR grid is lon centers -179.95..179.95 in 1..3600; lat centers 89.95..-89.95 */
        /* first pix represents -180..-179.9,90..89.9 */
        /* flip and add wrapping border to simplify interpolation */
        /* new grid is -180.1,180.1 in i=0,3601 and -90.1,90.1 in j=0,1801    */

        for (j = 0; j < ny; j++) {
            /* stmp points to ny x nx array read in from file */
            jl = j * (nx + 2); /* offset to line j in stmp */
            /* sref points to ny+2 * nx+2 array */
            /* flip so need to calculate line in sref */
            /* stmp line 0 goes into sref line ny */
            /* stmp line ny-1 goes into sref line 1 */
            jj = ny - j;
            jjl = jj * (nx + 2); /* offset to line jj in sref */
            for (i = 0; i < nx; i++) {
                if (stmp[jl + i] > -998.0) {
                    sref[jjl + i + 1] = stmp[jl + i];
                } else {
                    sref[jjl + i + 1] = sstbad;
                }
            }
            sref[jl] = sref[jl + nx];
            sref[jl + nx + 2 + nx + 1] = sref[jl + nx + 2 + 1];
        }
        ll = (ny + 1)*(nx + 2); /* start of last line */
        nl = ny * (nx + 2); /* start of next to last line */
        for (i = 0; i < nx + 2; i++) {
            sref[i] = sref[nx + 2 + i];
            sref[ll + i] = sref[nl + i];
        }
        free(stmp);
    }


    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), OI1NX + 1), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), OI1NY + 1), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    ftmp = 0.0;
    ntmp = 0;
    jl = (j - 1)*(nx + 2); /* offset to start of line j */
    nl = j * (nx + 2); /* offset to start of line j+1 */
    if (sref[jl + i ] > sstbad + 1) {
        ftmp += sref[jl + i ];
        ntmp++;
    }
    if (sref[jl + i + 1] > sstbad + 1) {
        ftmp += sref[jl + i + 1];
        ntmp++;
    }
    if (sref[nl + i + 1] > sstbad + 1) {
        ftmp += sref[nl + i + 1];
        ntmp++;
    }
    if (sref[nl + i ] > sstbad + 1) {
        ftmp += sref[nl + i ];
        ntmp++;
    }
    if (ntmp > 0) {
        ftmp /= ntmp;
        reftmp[0][0] = (sref[jl + i ] > sstbad + 1 ? sref[jl + i ] : ftmp);
        reftmp[0][1] = (sref[jl + i + 1] > sstbad + 1 ? sref[jl + i + 1] : ftmp);
        reftmp[1][1] = (sref[nl + i + 1] > sstbad + 1 ? sref[nl + i + 1] : ftmp);
        reftmp[1][0] = (sref[nl + i ] > sstbad + 1 ? sref[nl + i ] : ftmp);

        sst = (1 - t)*(1 - u) * reftmp[0][0] + t * (1 - u) * reftmp[0][1] + t * u * reftmp[1][1] + (1 - t) * u * reftmp[1][0];

    } else
        sst = sstbad;

    return (sst);
}

/* ------------------------------------------------------------------------------------ */
/* get_sstref() - retrieves reference sea surface temperature                           */
/*      ssttype has to be set for each pixel because climatology could be used          */
/*      instead of sstfile                                                              */
/*                                                                                      */
/* B. Franz, SAIC, May 2004.                                                           */
/* ---------------------------------------------------------------------------------------------- */
#include "smi_climatology.h"
#include "l1.h"

float get_sstref(short reftyp, char *sstfile, l1str *l1rec, int32_t ip) {
    float lon = l1rec->lon[ip];
    float lat = l1rec->lat[ip];
    float solz = l1rec->solz[ip];
    int32_t xsatid = l1rec->l1file->subsensorID;
    int32_t sensorID = l1rec->l1file->sensorID;
    int16_t yr, dy;
    double sec;

    unix2yds(l1rec->scantime, &yr, &dy, &sec);
    int32_t year = (int32_t) yr;
    int32_t day = (int32_t) dy;
    float sst = sstbad;
    int32_t sd_id;

    if (format < 0) {

        /* Does the file exist? */
        if (access(sstfile, F_OK) || access(sstfile, R_OK)) {
            printf("-E- %s line %d: SST input file '%s' does not exist or cannot open.\n",
                    __FILE__, __LINE__, sstfile);
            exit(1);
        }

        /* What is it? */
        if (reftyp == 8) {
            /* It's a 3-day Day or Night WindSat file from avgwindsat */
            format = WINDSAT3DN;

        } else if (reftyp == 7) {
            /* It's a 3-day WindSat file from SSMI */
            format = WINDSAT3DAY;

        } else if (reftyp == 6) {
            /* It's a daily WindSat file from SSMI */
            format = WINDSATDAY;

        } else if (reftyp == 5) {
            /* It's a 3-day Day or Night AMSR-E file from SSMI */
            format = AMSRE3DN;

        } else if (reftyp == 4) {
            /* It's a Daily NTEV2 file from Reynolds */
            format = NTEV2;

        } else if (reftyp == 3) {
            /* It's a Monthly ATSR file from Reynolds */
            format = ATSR;

        } else if (reftyp == 2) {
            /* It's a 3-day AMSR-E file from SSMI */
            format = AMSRE3DAY;

        } else if (reftyp == 1) {
            /* It's a daily AMSR-E file from SSMI */
            format = AMSREDAY;

        } else {

            char title[255] = "";

            // SDstart seg faults on the mac, so need to protect the 
            // call with Hishdf()
            if (Hishdf(sstfile))
                sd_id = SDstart(sstfile, DFACC_RDONLY);
            else
                sd_id = FAIL;

            if (sd_id != FAIL) {
                /* Format is HDF-like */
                char title[255] = "";
                if (SDreadattr(sd_id, SDfindattr(sd_id, "title"), (VOIDP) title) == 0) {
                    if (strstr(title, "Daily-OI-V2") != NULL) {
                        format = OISSTV2D;
                    } else if (strstr(title, "ARC Sea Surface Temperature") != NULL) {
                        format = ATSRDAY;

                    } else {
                        printf("-E- %s line %d: Unable to initialize SST file\n", __FILE__, __LINE__);
                        printf("%s %d\n", sstfile, day);
                        exit(1);
                    }
                } else if (SDreadattr(sd_id, SDfindattr(sd_id, "Title"), (VOIDP) title) == 0) {
                    if (strstr(title, "SST Climatology") != NULL) {
                        format = PATHCLIM;
                        if (smi_climatology_init(sstfile, day, SST) != 0) {
                            printf("-E- %s line %d: Unable to initialize SST file\n", __FILE__, __LINE__);
                            printf("%s %d\n", sstfile, day);
                            exit(1);
                        }
                    }
                }
                SDend(sd_id);

            } else {
                if (nc_open(sstfile, NC_NOWRITE, &sd_id) == NC_NOERR) {
                    /* Format is netCDF */
                    if (nc_get_att_text(sd_id, NC_GLOBAL, "title", title) == NC_NOERR) {
                        if (strstr(title, "Daily-OI") ||
                            strstr(title, "NOAA/NCEI 1/4 Degree Daily Optimum Interpolation Sea Surface Temperature (OISST) Analysis")) {
                            format = OISSTV2D;
                        } else if (strstr(title, "CMC") != NULL) {
                            format = CMC;
                        } else if (strstr(title, "Standard Mapped Image") != NULL) {
                            format = L3M;
                        }else {
                            printf("-E- %s : Unable to initialize SST file\n", __FILE__);
                            printf("%s %d\n", sstfile, day);
                            exit(1);
                        }
                    }
                    nc_close(sd_id);
                } else {

                    /* Format is not HDF or netCDF.  Assuming OISST Binary. */
                    format = OISSTBIN;
                }
            }
        }

        /* need a back-up climatology for NRT products */
        if (format != PATHCLIM) {
            char *tmp_str;
            char file [FILENAME_MAX] = "";
            if ((tmp_str = getenv("OCDATAROOT")) == NULL) {
                printf("OCDATAROOT environment variable is not defined.\n");
                exit(1);
            }
            strcpy(file, tmp_str);
            strcat(file, "/common/sst_climatology.hdf");
            if (smi_climatology_init(file, day, SST) != 0) {
                printf("-E- %s line %d: Unable to initialize SST file\n", __FILE__, __LINE__);
                printf("%s %d\n", file, day);
                exit(1);
            }
        }
    }

    switch (format) {
    case PATHCLIM:
        sst = smi_climatology(lon, lat, SST);
        l1rec->ssttype[ip] = 0; /* SST Climatology */
        /* don't know if this is skin or bulk, but sst doesn't use these so don't care */
        break;
    case OISSTBIN:
        sst = get_oisst(sstfile, lon, lat);
        l1rec->ssttype[ip] = 1; /* Reynolds Weekly 1 deg OI */
        /* All Reynolds reference fields are bulk, sst calculation needs skin */
        sst = sst - 0.17;
        break;
    case OISSTV2D:
        sst = get_oisstv2d(sstfile, lon, lat);
        l1rec->ssttype[ip] = 2; /* Reynolds Daily 0.25-deg OI */
        /* All Reynolds reference fields are bulk, sst calculation needs skin */
        sst = sst - 0.17;
        break;
    case AMSRE3DAY:
        sst = get_amsre(sstfile, lon, lat, solz, xsatid, sensorID);
        l1rec->ssttype[ip] = 3; /* AMSR-E 3-Day 0.25-deg */
        /* AMSR is microwave skin, needs correction to skin */
        sst = sst - 0.15;
        break;
    case AMSREDAY:
        sst = get_amsre(sstfile, lon, lat, solz, xsatid, sensorID);
        l1rec->ssttype[ip] = 4; /* AMSR-E Daily 0.25-deg */
        /* AMSR is microwave skin, needs correction to skin */
        sst = sst - 0.15;
        break;
    case ATSR:
        sst = get_atsr(sstfile, lon, lat, year, day, 1);
        l1rec->ssttype[ip] = 5; /* ATSR Monthly 0.25-deg */
        /* atsr is skin, no correction needed */
        break;
    case NTEV2:
        sst = get_ntev2(sstfile, lon, lat, year, day);
        l1rec->ssttype[ip] = 6; /* NTEV2 Daily 0.25-deg */
        /* All Reynolds reference fields are bulk, sst calculation needs skin */
        sst = sst - 0.17;
        break;
    case AMSRE3DN:
        sst = get_amsre(sstfile, lon, lat, solz, xsatid, sensorID);
        l1rec->ssttype[ip] = 7; /* AMSR-E 3-Day Day or Night 0.25-deg */
        /* AMSR is microwave skin, needs correction to skin */
        if (sst > sstbad + 1) {
            sst = sst - 0.15;
        }
        break;
    case ATSRDAY:
        sst = get_atsrday(sstfile, lon, lat);
        l1rec->ssttype[ip] = 8; /* ATSR Daily 0.1-deg */
        /* We use the sst_skin value from the ATSR files */
        /* convert Kelvin to C */
        sst = sst - 273.15;
        break;
    case WINDSAT3DAY:
        sst = get_windsat(sstfile, lon, lat, solz, xsatid, sensorID);
        l1rec->ssttype[ip] = 9; /* WindSat 3-Day 0.25-deg */
        /* WindSat is microwave skin, needs correction to skin */
        sst = sst - 0.15;
        break;
    case WINDSATDAY:
        sst = get_windsat(sstfile, lon, lat, solz, xsatid, sensorID);
        l1rec->ssttype[ip] = 10; /* WindSat Daily 0.25-deg */
        /* WindSat is microwave skin, needs correction to skin */
        sst = sst - 0.15;
        break;
    case WINDSAT3DN:
        sst = get_windsat(sstfile, lon, lat, solz, xsatid, sensorID);
        l1rec->ssttype[ip] = 11; /* WindSat 3-Day Day or Night 0.25-deg */
        /* WindSat is microwave skin, needs correction to skin */
        if (sst > sstbad + 1) {
            sst = sst - 0.15;
        }
        break;
    case CMC:
        sst = get_cmcsst(sstfile, lon, lat, 1); // units of Kelvin
        l1rec->ssttype[ip] = 12; /* CMC analyzed SST */
        /* CMC analyzed SST fields are "foundation" i.e. bulk, sst calculation needs skin correction*/
        sst = sst - 0.17;
        break;
    case L3M:
        sst = get_cmcsst(sstfile, lon, lat, 0);// units of degress_C
        l1rec->ssttype[ip] = 13; /* L3M  SST */
        break;
    default:
        printf("-E- %s line %d: unknown SST input file format for %s.\n", __FILE__, __LINE__, sstfile);
        break;
    }

    if (sst < (sstbad + 1)) {
        sst = smi_climatology(lon, lat, SST);
        l1rec->ssttype[ip] = 0; /* SST Climatology */
        /* don't know if this is skin or bulk, but sst doesn't use these so don't care */
    }

    return (sst);
}
