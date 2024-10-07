/* ============================================================================ */
/* module l1_l1c.c - functions to read L1C for l2gen (l1info really)            */
/* Written By:  Don Shea                                                        */
/*                                                                              */
/* ============================================================================ */

#include <netcdf.h>
#include "l1_l1c_anc.h"

#include <nc4utils.h>
#include "libnav.h"
#include <stdio.h>
#include <math.h>


static size_t num_scans, num_pixels;
static int ncid_L1C;

// time
static double start_time;
static double end_time;
static double delta_time;

// geolocation data
static int lonId, latId;
static float latFillValue = BAD_FLT;
static float lonFillValue = BAD_FLT;

/**
 * Open the OCI L1C Ancillary file and perform some one-time tasks (as opposed to tasks that are per scan)
 *
 * Get   
 * @param file
 * @return 
 */
int openl1_l1c_anc(filehandle * file) {
    int dimid, status;
    
    // Open the netcdf4 input file
    printf("Opening L1C file\n");
    status = nc_open(file->name, NC_NOWRITE, &ncid_L1C);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    
    // num_scans
    status = nc_inq_dimid(ncid_L1C, "bins_along_track", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading bins_along_track.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1C, dimid, &num_scans);

    // num_pixels
    status = nc_inq_dimid(ncid_L1C, "bins_across_track", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading bins_across_track.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1C, dimid, &num_pixels);

    if (want_verbose) {
        printf("L1C Ancillary Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);
    } // want_verbose

    // get start time
    size_t att_len;
    status = nc_inq_attlen(ncid_L1C, NC_GLOBAL, "time_coverage_start", &att_len);
    check_err(status, __LINE__, __FILE__);

    // allocate required space before retrieving values 
    char* time_str = (char *) malloc(att_len + 1); // + 1 for trailing null 

    // get attribute values 
    status = nc_get_att_text(ncid_L1C, NC_GLOBAL, "time_coverage_start", time_str);
    check_err(status, __LINE__, __FILE__);
    time_str[att_len] = '\0';
    start_time = isodate2unix(time_str);
    free(time_str);

    // get end time
    status = nc_inq_attlen(ncid_L1C, NC_GLOBAL, "time_coverage_end", &att_len);
    check_err(status, __LINE__, __FILE__);

    // allocate required space before retrieving values 
    time_str = (char *) malloc(att_len + 1); // + 1 for trailing null 

    // get attribute values 
    status = nc_get_att_text(ncid_L1C, NC_GLOBAL, "time_coverage_end", time_str);
    check_err(status, __LINE__, __FILE__);
    time_str[att_len] = '\0';
    end_time = isodate2unix(time_str);
    free(time_str);

    delta_time = (end_time - start_time) / num_scans;

    // Setup geofile pointers
    status = nc_inq_varid(ncid_L1C, "longitude", &lonId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(ncid_L1C, lonId, NULL, &lonFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(ncid_L1C, "latitude", &latId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(ncid_L1C, latId, NULL, &latFillValue);
    check_err(status, __LINE__, __FILE__);

    file->sd_id = ncid_L1C;
    file->npix = num_pixels;
    file->nscan = num_scans;
    file->ndets = 1;
    file->terrain_corrected = 1; // presumed.
    strcpy(file->spatialResolution, "5.2km");

    return (LIFE_IS_GOOD);
}


/**
 * Read the specified scan line from the specified L1B file. 
 * For each scan, get scan_delta_time_ms. 
 * Store Lt in l1rec->Lt.
 * Read GEO file.
 * 
 * 
 * @param file
 * @param line
 * @param l1rec
 * @return 
 */
int readl1_l1c_anc(filehandle *file, int32_t line, l1str *l1rec) {

    int status;
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    for (int ip = 0; ip < num_pixels; ip++) {
        l1rec->pixnum[ip] = ip;
    }

    l1rec->npix = file->npix;
    l1rec->scantime = start_time + (delta_time * line);

    int16_t syear, sday;
    double secs;
    unix2yds(l1rec->scantime, &syear, &sday, &secs);

    int32_t yr = syear;
    int32_t dy = sday;
    int32_t msec = (int32_t) (secs * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);

    l1rec->fsol = pow(1.0 / esdist, 2);

    start[0] = line;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = num_pixels;               // 1 line at a time
    count[2] = 1;
    status = nc_get_vara_float(ncid_L1C, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(ncid_L1C, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);

    for(int i=0; i<num_pixels; i++) {
        if(l1rec->lat[i] == latFillValue)
            l1rec->lat[i] = BAD_FLT;
        if(l1rec->lon[i] == lonFillValue)
            l1rec->lon[i] = BAD_FLT;

        // set solz to 0 so all pixels are daytime
        l1rec->solz[i] = 0.0;

    }
    
    return (LIFE_IS_GOOD);
}


/**
 * Close L1B file, GEO file, and free memory
 * @param file
 * @return 
 */
int closel1_l1c_anc(filehandle *file) {
    int status;

    printf("Closing L1C Ancillary file\n");
    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);

    return (LIFE_IS_GOOD);
}





