#include <sys/types.h>
#include <unistd.h>

/* ============================================================================ */
/* module sssref.c - retrieve salinity from WOA climatology                     */
/* ============================================================================ */

#include "l12_proto.h"
#include "h5io.h"

static float sssbad = BAD_FLT;
static int32_t format = -1;


/* ----------------------------------------------------------------------------------- */
/* get_woasssclim() - read and interpolate WOA salinity climatology                    */
/* ----------------------------------------------------------------------------------- */

#define WOASSSNX 360
#define WOASSSNY 180

float get_woasssclim(char *sssfile, float lon, float lat, int day) {
    static int firstCall = 1;
    static int nx = WOASSSNX;
    static int ny = WOASSSNY;
    static float dx = 360.0 / WOASSSNX;
    static float dy = 180.0 / WOASSSNY;

    typedef float sssref_t[WOASSSNX + 2];
    static sssref_t* sssref;
    static sssref_t* sss_sav;
    //static float sssref[WOASSSNY+2][WOASSSNX+2];
    //static float sss_sav[WOASSSNY+2][WOASSSNX+2];

    float sss = sssbad;
    int i, j, ii;
    float xx, yy;
    float t, u;

    if (firstCall) {
        sssref = (sssref_t*) malloc((WOASSSNY + 2) * sizeof (sssref_t));
        sss_sav = (sssref_t*) malloc((WOASSSNY + 2) * sizeof (sssref_t));

        float ssstmp[WOASSSNY][WOASSSNX];
        char name [H4_MAX_NC_NAME] = "";
        char sdsname[H4_MAX_NC_NAME] = "";
        int32 sd_id;
        int32 sds_id;
        int32 rank;
        int32 nt;
        int32 dims[H4_MAX_VAR_DIMS];
        int32 nattrs;
        int32 start[4] = {0, 0, 0, 0};
        int32 status, ix, iy, ct;
        int32 lymin, lxmin, lymax, lxmax;
        float sum;

        int16 mon = (int) day / 31 + 1; // month of year (no need for perfection)

        firstCall = 0;

        printf("Loading SSS reference from Climatology file: %s\n", sssfile);
        printf("\n");

        /* Open the file */
        sd_id = SDstart(sssfile, DFACC_RDONLY);
        if (sd_id == FAIL) {
            fprintf(stderr, "-E- %s line %d: error reading %s.\n",
                    __FILE__, __LINE__, sssfile);
            exit(1);
        }

        /* Get the SDS index */
        sprintf(sdsname, "sss%2.2i", mon);
        sds_id = SDselect(sd_id, SDnametoindex(sd_id, sdsname));
        status = SDgetinfo(sds_id, name, &rank, dims, &nt, &nattrs);
        status = SDreaddata(sds_id, start, NULL, dims, (VOIDP) & ssstmp[0][0]);
        if (status != 0) {
            fprintf(stderr, "-E- %s Line %d:  Error reading SDS %s from %s.\n",
                    __FILE__, __LINE__, sdsname, sssfile);
            exit(1);
        }
        status = SDendaccess(sds_id);
        status = SDend(sd_id);

        /* rotate 180-deg and add wrapping border to simplify interpolation */
        /* new grid is -180.5,180.5 in i=0,361 and -90.5,90.5 in j=0,181    */

        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                ii = (i < nx / 2) ? i + nx / 2 : i - nx / 2;
                if (ssstmp[j][i] > 0)
                    sssref[j + 1][ii + 1] = ssstmp[j][i];
                else
                    sssref[j + 1][ii + 1] = sssbad;
            }
            sssref[j + 1][0] = sssref[j + 1][nx];
            sssref[j + 1][nx + 1] = sssref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sssref[0] [i] = sssref[1][i];
            sssref[ny + 1][i] = sssref[ny][i];
        }
        /*
         *  Extend the grid by 1 point (use only full grid boxes later)
         */
        memcpy(sss_sav, sssref,
                (WOASSSNY + 2) * (WOASSSNX + 2) * sizeof ( float));
        for (iy = 0; iy < ny + 2; iy++) {
            lymin = (iy == 0) ? 0 : iy - 1;
            lymax = (iy == ny + 1) ? ny + 1 : iy + 1;

            for (ix = 0; ix < nx + 2; ix++) {
                if (sssref[iy][ix] < sssbad + 1) {
                    lxmin = (ix == 0) ? 0 : ix - 1;
                    lxmax = (ix == nx + 1) ? nx + 1 : ix + 1;

                    sum = 0.;
                    ct = 0;
                    for (j = lymin; j < lymax + 1; j++)
                        for (i = lxmin; i < lxmax + 1; i++) {
                            if (sss_sav[j][i] > sssbad + 1) {
                                sum += sss_sav[j][i];
                                ct++;
                            }
                        }
                    if (ct > 0)
                        sssref[iy][ix] = sum / ct;
                }
            }
        }
    } /* end 1-time grid set up */

    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), WOASSSNX), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), WOASSSNY), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, only using full grid boxes */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    if ((sssref[j][i] > sssbad + 1) && (sssref[j][i + 1] > sssbad + 1) &&
            (sssref[j + 1][i] > sssbad + 1) && (sssref[j + 1][i + 1] > sssbad + 1)) {

        sss = (1 - t)*(1 - u) * sssref[j ][i ]
                + t * (1 - u) * sssref[j ][i + 1]
                + t * u * sssref[j + 1][i + 1]
                + (1 - t) * u * sssref[j + 1][i ];

    } else
        sss = sssbad;

    return (sss);
}

/*
 *  for the HYCOM file usage...
 */
#define HYCOMNX 1440
#define HYCOMNY 720

float get_hycom_sss(char *sssfile, float lon, float lat, float *sss)
/*******************************************************************

   get_hycom_sss

   purpose: read in the daiy HYCOM sss and give the interpolated value 
     for the input lat, lon

   Returns type: int - return status: 0 is good, 1 if any other error
      will fail with -1 status if an h5 io error found (not associated with 
      the file not being h5 or correct format of h5)

   Parameters: (in calling order)
      Type              Name            I/O     Description
      ----              ----            ---     -----------
      char *            sssfile          I      sss file name
      float             lon              I      longitude: -180 -> 180
      float             lat              I      latitude: -90 -> 90
      float *           sss              O      sss found - will be sssbad 
                                                  if over land or otherwise 
                                                  bad value

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 12 Oct 2012     adapted from get_woasssclim in this file
      W. Robinson, SAIC 16 Oct 2012     this version 2 will interpolate using
                                        0 for any missing data

  Note that the HYCOM file only contains a dataset of salinity, so 
  the size of 1440, 720 probably should be an ID requirement as is the 
  existance of the dataset name 'salinity'

 *******************************************************************/ {
    static int firstCall = 1;
    static int nx = HYCOMNX;
    static int ny = HYCOMNY;
    static float dx = 360.0 / HYCOMNX;
    static float dy = 180.0 / HYCOMNY;

    typedef float sssref_t[ HYCOMNX + 2 ];
    static sssref_t* sssref;
    static sssref_t* sss_sav;
    //static float sssref[ HYCOMNY + 2 ][ HYCOMNX + 2 ];
    //static float sss_sav[ HYCOMNY + 2 ][ HYCOMNX + 2 ];

    int i, j, iter_ct, iter_max;
    float t, u, xx, yy;

    *sss = sssbad;
    if (firstCall) {
        sssref = (sssref_t*) malloc((HYCOMNY + 2) * sizeof (sssref_t));
        sss_sav = (sssref_t*) malloc((HYCOMNY + 2) * sizeof (sssref_t));
        float ssstmp[ny][nx], sum;
        h5io_str fid, dsid, grpid;
        int nobj, *types, ndim, dim_siz[10], sto_len, ix, iy, ct;
        int lymin, lxmin, lymax, lxmax;
        char **o_names;
        H5T_class_t class;
        hid_t native_typ;

        firstCall = 0;

        printf("Loading SSS reference from HYCOM file: %s\n", sssfile);
        printf("\n");

        /*
         *  check the file to see that it is correct - just return to caller if 
         *  not right format, but exit if something that should be performed can't be
         */
        if (h5io_openr(sssfile, 0, &fid) != 0)
            return 1;
        /* the group may need setting to work, see if '/' will do */
        if (h5io_set_grp(&fid, "/", &grpid) != 0) {
            fprintf(stderr, "-E- %s, %d: Failed to set trivial group\n",
                    __FILE__, __LINE__);
            exit(-1);
        }
        if (h5io_grp_contents(&grpid, &nobj, &o_names, &types) != 0) {
            fprintf(stderr, "-E- %s, %d: Failed to find contents of sss file\n",
                    __FILE__, __LINE__);
            exit(-1);
        }
        /* first name should be salinity */
        if (strcmp(o_names[0], "salinity") != 0)
            return 1;

        /* set to the salinity dataset and get info */
        if (h5io_set_ds(&fid, "salinity", &dsid) != 0)
            return 1;
        if (h5io_info(&dsid, NULL, &class, &native_typ, &ndim, dim_siz,
                &sto_len) != 0) {
            fprintf(stderr, "-E- %s, %d: Failed to find info on salinity dataset\n",
                    __FILE__, __LINE__);
            exit(-1);
        }

        if ((ndim != 2) || (dim_siz[0] != 720) || (dim_siz[1] != 1440))
            return 1;
        /*
         *  at this point, we have identified the salinity file, just read 
         *  the dataset now
         */
        if (h5io_rd_ds(&dsid, (void *) ssstmp) != 0) {
            fprintf(stderr, "-E- %s, %d: Failed to read salinity dataset\n",
                    __FILE__, __LINE__);
            exit(-1);
        }
        /* close and go */
        h5io_close(&dsid);
        h5io_close(&fid);

        /* add wrapping border to simplify interpolation */
        /* new grid is -180.125,180.125 in i=0,1441 and -90.125,90.125 in j=0,721 */
        for (j = 0; j < ny; j++) {
            for (i = 0; i < nx; i++) {
                if ((ssstmp[j][i] > 0) && (ssstmp[j][i] <= 1000.))
                    sssref[j + 1][i + 1] = ssstmp[j][i];
                else
                    sssref[j + 1][i + 1] = sssbad;
            }
            sssref[j + 1][0] = sssref[j + 1][nx];
            sssref[j + 1][nx + 1] = sssref[j + 1][1];
        }
        for (i = 0; i < nx + 2; i++) {
            sssref[0] [i] = sssref[1][i];
            sssref[ny + 1][i] = sssref[ny][i];
        }
        /*
         *  extend data to have good interpolation points at edge = coast
         */
        iter_ct = 0;
        iter_max = 4;
        do {
            iter_ct++;
            memcpy(sss_sav, sssref,
                    (HYCOMNY + 2) * (HYCOMNX + 2) * sizeof ( float));
            for (iy = 0; iy < ny + 2; iy++) {
                lymin = (iy == 0) ? 0 : iy - 1;
                lymax = (iy == ny + 1) ? ny + 1 : iy + 1;

                for (ix = 0; ix < nx + 2; ix++) {
                    if (sssref[iy][ix] < sssbad + 1) {
                        lxmin = (ix == 0) ? 0 : ix - 1;
                        lxmax = (ix == nx + 1) ? nx + 1 : ix + 1;

                        sum = 0.;
                        ct = 0;
                        for (j = lymin; j < lymax + 1; j++)
                            for (i = lxmin; i < lxmax + 1; i++) {
                                if (sss_sav[j][i] > sssbad + 1) {
                                    sum += sss_sav[j][i];
                                    ct++;
                                }
                            }
                        if (ct > 0)
                            sssref[iy][ix] = sum / ct;
                    }
                }
            }
        } while (iter_ct < iter_max);
        /*  and end set up */
    }

    /* locate LL position within reference grid */
    i = MAX(MIN((int) ((lon + 180.0 + dx / 2) / dx), nx), 0);
    j = MAX(MIN((int) ((lat + 90.0 + dy / 2) / dy), ny), 0);

    /* compute longitude and latitude of that grid element */
    xx = i * dx - 180.0 - dx / 2;
    yy = j * dy - 90.0 - dy / 2;

    /* bilinearly interpolate, replacing missing (land) values with 
        average of valid values in box */
    t = (lon - xx) / dx;
    u = (lat - yy) / dy;

    if ((sssref[j][i] > sssbad + 1) && (sssref[j][i + 1] > sssbad + 1) &&
            (sssref[j + 1][i] > sssbad + 1) && (sssref[j + 1][i + 1] > sssbad + 1))
        *sss = (1 - t)*(1 - u) * sssref[j ][i ]
            + t * (1 - u) * sssref[j ][i + 1]
            + t * u * sssref[j + 1][i + 1]
            + (1 - t) * u * sssref[j + 1][i ];
    else
        *sss = sssbad;

    return 0;
}

/* ----------------------------------------------------------------------------------- */
/* get_sssref() - retrieves reference sea surface salinity                 .           */
/* ----------------------------------------------------------------------------------- */
#define WOACLIM 1

float get_sssref(char *sssfile, float lon, float lat, int day) {
    float sss = sssbad;
    int32_t sd_id;

    if (format < 0) {

        /* Does the file exist? */
        if (access(sssfile, F_OK) || access(sssfile, R_OK)) {
            printf("-E- %s: SSS input file '%s' does not exist or cannot open.\n",
                    __FILE__, sssfile);
            exit(1);
        }

        /* What is it? */
        sd_id = SDstart(sssfile, DFACC_RDONLY);
        if (sd_id != FAIL) {
            /* Format is HDF-like */
            char title[255] = "";
            if (SDreadattr(sd_id, SDfindattr(sd_id, "Title"), (VOIDP) title) == 0) {
                if (strstr(title, "WOA Sea Surface Salinity Monthly Climatology") != NULL) {
                    format = WOACLIM;
                }
            } else {
                printf("-E- %s: unable to read SSS input file %s.\n", __FILE__, sssfile);
                exit(1);
            }
            SDend(sd_id);
        }
        else {
            if (get_hycom_sss(sssfile, lon, lat, &sss) == 0)
                format = 2;
        }
    }

    switch (format) {
    case WOACLIM:
        sss = get_woasssclim(sssfile, lon, lat, day);
        break;
    case 2:
        get_hycom_sss(sssfile, lon, lat, &sss);
        break;
    default:
        printf("-E- %s: unknown SSS input file format for %s.\n", __FILE__, sssfile);
        exit(1);
        break;
    }

    // assume fresh water if no data

    if (sss < 0.0) {
        sss = 0.0;
    }

    return (sss);
}


