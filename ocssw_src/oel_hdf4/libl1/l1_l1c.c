/* ============================================================================ */
/* module l1_l1c.c - functions to read L1C for l2gen (l1info really)            */
/* Written By:  Don Shea                                                        */
/*                                                                              */
/* ============================================================================ */

#include <netcdf.h>
#include "l1_l1c.h"

#include <nc4utils.h>
#include "libnav.h"
#include <stdio.h>
#include <math.h>


static short *tmpShort;

static size_t num_scans, num_pixels;
static size_t num_bands;
static int ncid_L1C;

// time
static double *scan_time;
static double file_start_day;

// geolocation data
static int geolocationGrp;
static int lonId, latId;
static float latFillValue = BAD_FLT;
static float lonFillValue = BAD_FLT;

/**
 * Open the OCI L1B file and perform some one-time tasks (as opposed to tasks that
 * are per scan), including:
 *   -Get L1B dimensions num_scans, num_bands, num_pixels (static).
 *   locate memory for some static arrays based upon these dimensions.
 *   -Get L1B group ids e.g. l1bScanLineGrp, observationGrp and it's 8 "band_%d" var ids (static)
 *

 * Get   
 * @param file
 * @return 
 */
int openl1_l1c(filehandle * file) {
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

    // num_bands
    status = nc_inq_dimid(ncid_L1C, "intensity_bands_per_view", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading intensity_bands_per_view.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1C, dimid, &num_bands);

    if (want_verbose) {
        printf("L1C Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);
    } // want_verbose

    // allocate all of the data
    tmpShort = (short*) malloc(num_pixels * sizeof(short));
    scan_time = (double*) malloc(num_scans * sizeof(double));

    // Get group id from L1B file for GROUP bin_attributes.
    int groupid;
    if ((nc_inq_grp_ncid(ncid_L1C, "bin_attributes", &groupid)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding bin_attributes.\n");
        exit(EXIT_FAILURE);
    }  
    int varid;
    double scan_timeFillValue = BAD_FLT;
    status = nc_inq_varid(groupid, "nadir_view_time", &varid);
    if(status != NC_NOERR) {
        fprintf(stderr, "-E- Error finding nadir_view_time.\n");
        exit(EXIT_FAILURE);
    }  
    status = nc_inq_var_fill(groupid, varid, NULL, &scan_timeFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_double(groupid, varid, scan_time);
    check_err(status, __LINE__, __FILE__);    

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

    double start_time = isodate2unix(time_str);
    int16_t syear, smon, sday;
    double secs;
    unix2ymds(start_time, &syear, &smon, &sday, &secs);
    file_start_day = ymds2unix(syear, smon, sday, 0.0);

    free(time_str);

    for(int i=0; i<num_scans; i++) {
        if(scan_time[i] == scan_timeFillValue)
            scan_time[i] = BAD_FLT;
    }
    
    // Setup geofile pointers
    status = nc_inq_grp_ncid(ncid_L1C, "geolocation_data", &geolocationGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(geolocationGrp, "longitude", &lonId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, lonId, NULL, &lonFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(geolocationGrp, "latitude", &latId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, latId, NULL, &latFillValue);
    check_err(status, __LINE__, __FILE__);

    file->sd_id = ncid_L1C;
    file->nbands = num_bands;
    file->npix = num_pixels;
    file->nscan = num_scans;
    file->ndets = 1;
    file->terrain_corrected = 1; // presumed.
    strcpy(file->spatialResolution, "5.2km");

    if (want_verbose)
        printf("file->nbands = %d\n", (int) file->nbands);

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
int readl1_l1c(filehandle *file, int32_t line, l1str *l1rec) {

    int status;
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    for (int ip = 0; ip < num_pixels; ip++) {
        l1rec->pixnum[ip] = ip;
    }

    l1rec->npix = file->npix;
    if(scan_time[line] == BAD_FLT) {
        l1rec->scantime = BAD_FLT;
        l1rec->fsol = BAD_FLT;
    } else {
        l1rec->scantime = file_start_day + scan_time[line];

        int16_t syear, sday;
        double secs;
        unix2yds(l1rec->scantime, &syear, &sday, &secs);

        int32_t yr = syear;
        int32_t dy = sday;
        int32_t msec = (int32_t) (secs * 1000.0);
        double esdist = esdist_(&yr, &dy, &msec);

        l1rec->fsol = pow(1.0 / esdist, 2);
    }

    start[0] = line;
    start[1] = 0;
    start[2] = 0;
    count[0] = 1;
    count[1] = num_pixels;               // 1 line at a time
    count[2] = 1;
    status = nc_get_vara_float(geolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start, count, l1rec->lon);
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
int closel1_l1c(filehandle *file) {
    int status;

    printf("Closing L1C file\n");
    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);

    // Free memory
    // From openl1_oci
    if (tmpShort) free(tmpShort);
    if (scan_time) free(scan_time);

    return (LIFE_IS_GOOD);
}





