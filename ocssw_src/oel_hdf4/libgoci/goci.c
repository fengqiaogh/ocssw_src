/**
 *  @file libgoci/goci.c
 *  @brief GOCI file format reader
 *  @author Paul Martinolich
 *  @author Naval Research Laboratory
 */

#include <time.h>
#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <inttypes.h>
#include <unistd.h>
#include <proj.h>
#include <math.h>

#include "hdf5.h"
#include "hdf5_hl.h"
#include "goci.h"

/**
 *  @brief init proj4 for GOCI geolocation
 *  @param[in/out]  goci_l1b   GOCI file structure
 *
 *  Uses the information in the goci_l1b structure to init the proj4 structure
 *  for converting the geolocation information.
 *
 *  @returns -1 for open error, -2 for close error, -4 for read error
 */
int
goci_proj4_open(goci_l1b_t *l1b) {
    hid_t gid;
    hid_t attr;
    herr_t herr;
    float centralLat;
    float centralLon;
    float equitorialRadius;
    float polarRadius;
    float llLat, llLon;
    float urLat, urLon;
    char projStr[1024];

    // read Map Projection attributes
    gid = H5Gopen(l1b->fid, "/HDFEOS/POINTS/Map Projection", H5P_DEFAULT);
    if (gid < 0) {
        return -1;
    }
    attr = H5Aopen_name(gid, "Central Latitude (parallel)");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &centralLat);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    attr = H5Aopen_name(gid, "Central Longitude (meridian)");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &centralLon);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    attr = H5Aopen_name(gid, "Equitorial radius of Earth ellipsoid");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &equitorialRadius);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    attr = H5Aopen_name(gid, "Polar radius of Earth ellipsoid");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &polarRadius);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    herr = H5Gclose(gid);
    if (herr < 0)
        return -2;


    // read Scene Header attributes
    gid = H5Gopen(l1b->fid, "/HDFEOS/POINTS/Scene Header", H5P_DEFAULT);
    if (gid < 0) {
        return -1;
    }
    attr = H5Aopen_name(gid, "Scene lower-left latitude");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &llLat);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    attr = H5Aopen_name(gid, "Scene lower-left longitude");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &llLon);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    attr = H5Aopen_name(gid, "Scene upper-right latitude");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &urLat);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    attr = H5Aopen_name(gid, "Scene upper-right longitude");
    if (attr < 0) {
        H5Gclose(gid);
        return -1;
    }
    herr = H5Aread(attr, H5T_NATIVE_FLOAT, &urLon);
    H5Aclose(attr);
    if (herr < 0) {
        H5Gclose(gid);
        return -4;
    }

    PJ *pj;
    PJ_COORD c, c_ll, c_ur;
    
    // init the proj4 projections
    sprintf(projStr, "+proj=ortho +ellps=WGS84 +datum=WGS84 +lon_0=%g +lat_0=%g +a=%g +b=%g",
            centralLon, centralLat, equitorialRadius, polarRadius);
    pj = proj_create_crs_to_crs(PJ_DEFAULT_CTX,
                               projStr,
                               "+proj=longlat +ellps=WGS84 +datum=WGS84",
                               NULL);
    if(pj == NULL) {
        printf("Error - GOCI first PROJ projection failed to init\n");
        exit(1);
    }
    l1b->pj = proj_normalize_for_visualization(PJ_DEFAULT_CTX, pj);
    if(l1b->pj == NULL) {
        printf("Error - GOCI visualization PROJ projection failed to init\n");
        exit(1);
    }
    proj_destroy(pj);
    
    // calculate start and delta for grid
    // set default z and t
    c.xyzt.z = 0.0;
    c.xyzt.t = HUGE_VAL;

    c.xy.x = llLon;
    c.xy.y = llLat;
    c_ll = proj_trans(l1b->pj, PJ_INV, c);

    c.xy.x = urLon;
    c.xy.y = urLat;
    c_ur = proj_trans(l1b->pj, PJ_INV, c);

    l1b->startX = c_ll.xy.x;
    l1b->startY = c_ur.xy.y;
    l1b->deltaX = (c_ur.xy.x - c_ll.xy.x) / l1b->npixels;
    l1b->deltaY = (c_ll.xy.y - c_ur.xy.y) / l1b->nscans;

    return 0;

}

/**
 *  @brief free all memory related to proj4
 *  @param[in/out]  goci_l1b   GOCI file structure
 *
 *  Free up the memory used by the proj4 geolocation functions
 *
 */
void
goci_proj4_close(goci_l1b_t *l1b) {
    proj_destroy(l1b->pj);
}

/**
 *  @brief open GOCI bands
 *  @param[in]   src_path   path to one GOCI file
 *  @param[out]  goci_l1b   GOCI file structure
 *
 *  Opens a GOCI data file and attaches to the main data set "Image Data"
 *  which is stored as a GRID.
 *
 *  @returns -1 for open error, -2 for close error, -3 for invalid format
 *           -4 for read error, -5 for memory error
 *  W. Robinson, SAIC, Nov 2014, the position in the file is really ( radius, 
 *    sub-sat lon, sub-sat lat ), so derive position accordingly
 *  W. Robinson, SAIC, 5 Dec 2014, set up the slot navigation to get better
 *    times
 */

int
goci_l1b_open(const char *src_path, goci_l1b_t **goci_l1b) {
    goci_l1b_t *l1b = NULL;

    hsize_t dims[3];
    herr_t herr;
    htri_t exists;
    hid_t fid, gid;
    hid_t loc_id, dset_id[8];
    hid_t space, dset;
    float cpos[3], radius_in_xy;


    char *fields[] = {"Band 1 Image Pixel Values",
        "Band 2 Image Pixel Values",
        "Band 3 Image Pixel Values",
        "Band 4 Image Pixel Values",
        "Band 5 Image Pixel Values",
        "Band 6 Image Pixel Values",
        "Band 7 Image Pixel Values",
        "Band 8 Image Pixel Values"};
    int i;
    int32_t slot_nav_avail; // is the slot navigation good?
    static float slot_rel_time[16]; // relative time of each slot to start
    unsigned char *slot_asg; // array of slot numbers for scene pixels

    // open GOCI file

    fid = H5Fopen(src_path, H5F_ACC_RDONLY, H5P_DEFAULT);
    if (fid < 0)
        return -1;

    gid = H5Gopen(fid, "/", H5P_DEFAULT);
    if (gid < 0) {
        herr = H5Fclose(fid);
        if (herr < 0)
            return -2;
        return -1;
    }

    // open group that holds GRIDS fields;
    loc_id = -1;
    exists = H5Lexists(gid, "HDFEOS/GRIDS/Image Data/Data Fields", H5P_DEFAULT);
    if (exists > 0)
        loc_id = H5Gopen(gid, "HDFEOS/GRIDS/Image Data/Data Fields", H5P_DEFAULT);

    if (loc_id < 0)
        goto grid_error;

    // open each data set in the group
    for (i = 0; i < 8; i++) {

        dset = H5Dopen(loc_id, fields[i], H5P_DEFAULT);
        if (dset < 0)
            goto grid_error;

        dset_id[i] = dset;
    }
    H5Gclose(loc_id);

    // make object

    l1b = (goci_l1b_t *) malloc(sizeof (goci_l1b_t));
    if (l1b == NULL) {
        H5Gclose(gid);
        H5Fclose(fid);
        return -5;
    }

    l1b->fid = fid;
    l1b->gid = gid;
    //    l1b->loc_id = loc_id;
    for (i = 0; i < 8; i++)
        l1b->dset_id[i] = dset_id[i];

    // set file extents from first band

    space = H5Dget_space(dset_id[0]);
    if (space < 0)
        goto grid_error;
    H5Sget_simple_extent_dims(space, dims, NULL);

    l1b->nbands = 8;
    l1b->npixels = dims[1];
    l1b->nscans = dims[0];

    l1b->pj = NULL;

    // get satellite position at scene center
    herr = H5LTget_attribute_float(l1b->gid, "HDFEOS/POINTS/Ephemeris",
            "Satellite position XYZ (ECEF) at scene center time", cpos);
    if (herr < 0)
        goto grid_error;

    // convert ( radius, sub-sat lon, sub-sat lat ) to position
    // (note that lat, lon already in radians)
    l1b->sat_pos[2] = cpos[0] * sin(cpos[2]) / 1000.; // z
    radius_in_xy = cpos[0] * cos(cpos[2]) / 1000.;
    l1b->sat_pos[1] = radius_in_xy * sin(cpos[1]); // y
    l1b->sat_pos[0] = radius_in_xy * cos(cpos[1]); // x

    // set up goci slot navigation for time derivation
    if ((slot_asg = (unsigned char *)
            malloc(dims[0] * dims[1] * sizeof ( unsigned char))) == NULL) {
        printf("%s,%d:E Unable to allocate space for slot_asg array\n",
                __FILE__, __LINE__);
        goto grid_error;
    }
    if (goci_slot_init(fid, dims, slot_rel_time, slot_asg, &slot_nav_avail)
            != 0)
        goto grid_error;
    l1b->slot_nav_avail = slot_nav_avail;
    l1b->slot_rel_time = slot_rel_time;
    l1b->slot_asg = slot_asg;

    *goci_l1b = l1b;

    return 0;

    // We had an error reading the grid projection information

grid_error:
    if (l1b)
        free(l1b);
    herr = H5Gclose(gid);
    if (herr < 0) {
        H5Fclose(fid);
        return -2;
    }
    herr = H5Fclose(fid);
    if (herr < 0)
        return -2;
    return -3;
}

/**
 *  @brief close GOCI bands
 *  @param[in]   goci_l1b  GOCI file
 *
 *  Close the GOCI file
 *
 *  @returns 0 for success, -2 for close error
 */

int
goci_l1b_close(goci_l1b_t *goci_l1b) {
    herr_t herr;
    int i;

    // close all data sets

    for (i = 0; i < 8; i++)
        herr = H5Dclose(goci_l1b->dset_id[i]);
    if (herr < 0)
        return -2;

    // close GRIDS Fields group

    //    herr = H5Gclose( goci_l1b->loc_id );
    //    if ( herr < 0 )
    //        return -2;

    // close file

    herr = H5Gclose(goci_l1b->gid);
    if (herr < 0)
        return -2;
    herr = H5Fclose(goci_l1b->fid);
    if (herr < 0)
        return -2;

    free(goci_l1b);

    return 0;
}

/**
 *  @brief read a date of GOCI data
 *  @param[in]  goci_l1b    GOCI l1b file information structure
 *  @param[in]  tstr_name   name of time string to read from (char *)
 *  @param[out] year        Year
 *  @param[out] month       Month (1-12)
 *  @param[out] day         Day (1-31)
 *
 *  This function gets the date from the GOCI data file.
 *
 *  @returns 0 for success, -3 for invalid format
 */

int
goci_l1b_get_date(goci_l1b_t *goci_l1b, char *tstr_name, int *year,
        int *month, int *day) {

    struct tm time;
    herr_t herr;

    char str[32];

    // For example, 12-APR-2011 02:15:38.398
    herr = H5LTget_attribute_string(goci_l1b->gid, "HDFEOS/POINTS/Ephemeris",
            tstr_name, str);
    if (herr < 0)
        return -3;

    strptime(str, "%d-%B-%Y %H:%M:%S.%f", &time);

    *year = time.tm_year + 1900;
    *month = time.tm_mon + 1;
    *day = time.tm_mday;

    return 0;
}

/**
 *  @brief read a time of GOCI data
 *  @param[in]  goci_l1b    GOCI file
 *  @param[in]  tstr_name   name of time string to read from (char *)
 *  @param[out] hour        Hour of day (0-23)
 *  @param[out] min         Minute (0-59)
 *  @param[out] sec         Second (0-59)
 *
 *  This function gets the time from the GOCI data file.
 *
 *  @returns 0 for success, -3 for invalid format
 */

int
goci_l1b_get_time(goci_l1b_t *goci_l1b, char *tstr_name, int *hour, int *min,
        int *sec) {
    struct tm time;
    herr_t herr;

    char str[32];

    // For example, 12-APR-2011 02:15:38.398
    herr = H5LTget_attribute_string(goci_l1b->gid, "HDFEOS/POINTS/Ephemeris", tstr_name, str);
    if (herr < 0) {
        fprintf(stderr, "H5Lt failed\n");
        return -3;
    }
    strptime(str, "%d-%B-%Y %H:%M:%S.%f", &time);

    *hour = time.tm_hour;
    *min = time.tm_min;
    *sec = time.tm_sec;

    return 0;
}

/**
 *  @brief read specified band and line from GOCI file.
 *  @param[in]  goci_l1b    GOCI file
 *  @param[in]  band        Band number to read (0-7)
 *  @param[in]  line        Line to read
 *  @param[out] buf         buffer to contain data
 *
 *  This function reads a scan line from the GOCI data.  The returned
 *  value can be converted to radiance by multiply the value by 1E-07.
 *
 *  @returns 0 for success, -4 for read error, -6 for invalid parameter
 */

int
goci_l1b_get_band(goci_l1b_t *goci_l1b, int band, int line, uint32_t *buf) {
    hsize_t start[2], count[2], dims[2];
    herr_t herr;
    hid_t space_file, space_mem;

    if (band < 0 || band > 7)
        return -6;

    // get data space of file
    space_file = H5Dget_space(goci_l1b->dset_id[band]);
    if (space_file < 0)
        return -4;
    H5Sget_simple_extent_dims(space_file, dims, NULL);

    start[0] = line;
    start[1] = 0;

    count[0] = 1;
    count[1] = dims[1];

    // define a memory space

    space_mem = H5Screate_simple(2, count, NULL);
    if (space_mem < 0)
        return -4;
    H5Sselect_all(space_mem);

    // define hyperslab

    herr = H5Sselect_hyperslab(space_file, H5S_SELECT_SET, start, NULL, count, NULL);
    if (herr < 0) {
        fprintf(stderr, "error reading %d\n", line);
        return -4;
    }

    // read line

    herr = H5Dread(goci_l1b->dset_id[band], H5T_STD_U32LE, space_mem, space_file, H5P_DEFAULT, buf);
    if (herr < 0) {
        fprintf(stderr, "error reading %d\n", line);
        return -4;
    }

    H5Sclose(space_mem);
    return 0;
}
