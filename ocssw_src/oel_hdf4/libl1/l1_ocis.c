/* ============================================================================ */
/* module l1b_ocis.c - functions to read OCI Similated  L1B for l2gen           */
/* Written By:  Don Shea                                                        */
/*                                                                              */
/* ============================================================================ */

/* Issues:
        -Assume num_pixels = ccd_pixels = SWIR_pixels
*/

#include <netcdf.h>
#include "l1_ocis.h"

#include "l1.h"
#include <nc4utils.h>
#include "libnav.h"
#include <stdio.h>
#include <math.h>
#include <allocate2d.h>


static float *Fobar; // reflectance to radiance conversion factors
static int extract_pixel_start = 0;
static int normalizedLt = 1;
static int use_rhot = 0;

static short *tmpShort;

// whole file stuff
static size_t expected_num_blue_bands = 120;
static size_t expected_num_red_bands = 120;
static size_t expected_num_SWIR_bands = 9;

static size_t num_scans, num_pixels;
static size_t num_blue_bands, num_red_bands, num_SWIR_bands;
static size_t tot_num_bands = 239;
static int ncid_L1B;

// scan line attributes
static double *scan_time; // seconds of day
static int32_t scan_time_year, scan_time_month, scan_time_day;

//static uint_8 *scan_quality;

// geolocation data
static int geolocationGrp;
static int lonId, latId, senzId, senaId, solzId, solaId;
static float latFillValue = -999.0;
static float lonFillValue = -999.0;
static short senzFillValue = -32768;
static float senzScale = 0.01;
static float senzOffset = 0.0;
static short senaFillValue = -32768;
static float senaScale = 0.01;
static float senaOffset = 0.0;
static short solzFillValue = -32768;
static float solzScale = 0.01;
static float solzOffset = 0.0;
static short solaFillValue = -32768;
static float solaScale = 0.01;
static float solaOffset = 0.0;

// Observation data
static int observationGrp;
static int Lt_blueId, Lt_redId, Lt_SWIRId; 
static float **Lt_blue;                //[num_blue_bands][num_pixels], This scan
static float **Lt_red;                 //[num_red_bands][num_pixels], This scan
static float **Lt_SWIR;                //[num_SWIR_bands][num_pixels], This scan
static float Lt_blueFillValue = -999.0;
static float Lt_redFillValue = -999.0;
static float Lt_SWIRFillValue = -999.0;

/**
 * Open the OCI L1B file and perform some one-time tasks (as opposed to tasks that
 * are per scan), including:
 *   -Get L1B dimensions num_scans, num_bands, num_pixels (static).
 *    Allocate memory for some static arrays based upon these dimensions.
 *   -Get L1B group ids e.g. l1bScanLineGrp, observationGrp and it's 8 "band_%d" var ids (static)
 *

 * Get   
 * @param file
 * @return 
 */
int openl1_ocis(filehandle * file) {
    int dimid, status;
    
    // Open the netcdf4 input file
    printf("Opening OCIS L1B file\n");
    status = nc_open(file->name, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    
    // num_scans
    status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading number_of_scans.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

    // num_pixels
    status = nc_inq_dimid(ncid_L1B, "ccd_pixels", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_pixels.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);

    // num_blue_bands
    status = nc_inq_dimid(ncid_L1B, "blue_bands", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_blue_bands.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1B, dimid, &num_blue_bands);
    if(num_blue_bands < expected_num_blue_bands) {
        fprintf(stderr, "-E- Not enough blue bands, expecting %d, found %d.\n",
                (int)expected_num_blue_bands, (int)num_blue_bands);
        exit(EXIT_FAILURE);
    }

    // num_red_bands
    status = nc_inq_dimid(ncid_L1B, "red_bands", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_red_bands.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1B, dimid, &num_red_bands);
    if(num_red_bands < expected_num_red_bands) {
        fprintf(stderr, "-E- Not enough red bands, expecting %d, found %d.\n",
                (int)expected_num_red_bands, (int)num_red_bands);
        exit(EXIT_FAILURE);
    }

    // num_SWIR_bands
    status = nc_inq_dimid(ncid_L1B, "SWIR_bands", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_SWIR_bands.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1B, dimid, &num_SWIR_bands);
    if(num_SWIR_bands < expected_num_SWIR_bands) {
        fprintf(stderr, "-E- Not enough SWIR bands, expecting %d, found %d.\n",
                (int)expected_num_SWIR_bands, (int)num_SWIR_bands);
        exit(EXIT_FAILURE);
    }

    if (want_verbose) {
        printf("OCIS L1B Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);
    } // want_verbose

    // allocate all of the data
    tmpShort = (short*) malloc(num_pixels * sizeof(short));
    scan_time = (double*) malloc(num_scans * sizeof(double));

    Lt_blue = allocate2d_float(num_blue_bands, num_pixels);
    Lt_red = allocate2d_float(num_red_bands, num_pixels);
    Lt_SWIR = allocate2d_float(num_SWIR_bands, num_pixels);

    // Get group id from L1B file for GROUP scan_line_attributes.
    int groupid;
    if ((nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &groupid)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding scan_line_attributes.\n");
        exit(EXIT_FAILURE);
    }   
    int varid;
    double scan_timeFillValue = -999.9;
    status = nc_inq_varid(groupid, "time", &varid);
    if(status == NC_NOERR) {
        status = nc_inq_var_fill(groupid, varid, NULL, &scan_timeFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(groupid, varid, scan_time);
        check_err(status, __LINE__, __FILE__);    
        status = nc_get_att_int(groupid, varid, "year", &scan_time_year);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_att_int(groupid, varid, "month", &scan_time_month);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_att_int(groupid, varid, "day", &scan_time_day);
        check_err(status, __LINE__, __FILE__);
    } else {
        status = nc_inq_varid(groupid, "ev_mid_time", &varid);
        check_err(status, __LINE__, __FILE__);    
        status = nc_inq_var_fill(groupid, varid, NULL, &scan_timeFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(groupid, varid, scan_time);
        check_err(status, __LINE__, __FILE__);    

        // get start time
        size_t att_len;
        status = nc_inq_attlen(ncid_L1B, NC_GLOBAL, "time_coverage_start", &att_len);
        check_err(status, __LINE__, __FILE__);

        // allocate required space before retrieving values 
        char* time_str = (char *) malloc(att_len + 1); // + 1 for trailing null 

        // get attribute values 
        status = nc_get_att_text(ncid_L1B, NC_GLOBAL, "time_coverage_start", time_str);
        check_err(status, __LINE__, __FILE__);
        time_str[att_len] = '\0';

        double start_time = isodate2unix(time_str);
        int16_t syear, smon, sday;
        double secs;
        unix2ymds(start_time, &syear, &smon, &sday, &secs);
        scan_time_year = syear;
        scan_time_month = smon;
        scan_time_day = sday;

    }

    for(int i=0; i<num_scans; i++) {
        if(scan_time[i] == scan_timeFillValue)
            scan_time[i] = BAD_FLT;
    }
    
    // check if the normalizedLt attribute is set; GMAO output will be normalized Lt, and likley won't have 
    // this attribute set, so default to TRUE.  If the file defines it, use it.
    int normLt;
    status = nc_get_att_int(ncid_L1B, NC_GLOBAL, "normalizedLt", &normLt);
    if (!status) {
        normalizedLt = normLt;
    }

    // read the orbit#
    int orbit_number;
    status = nc_get_att_int(ncid_L1B, NC_GLOBAL, "orbit_number", &orbit_number);
    check_err(status, __LINE__, __FILE__);

    // Setup geofile pointers
    status = nc_inq_grp_ncid(ncid_L1B, "geolocation_data", &geolocationGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(geolocationGrp, "longitude", &lonId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, lonId, NULL, &lonFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(geolocationGrp, "latitude", &latId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, latId, NULL, &latFillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "sensor_zenith", &senzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, senzId, NULL, &senzFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senzId, "scale_factor", &senzScale);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senzId, "add_offset", &senzOffset);
    check_err(status, __LINE__, __FILE__);
    
    status = nc_inq_varid(geolocationGrp, "sensor_azimuth", &senaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, senaId, NULL, &senaFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senaId, "scale_factor", &senaScale);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senaId, "add_offset", &senaOffset);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "solar_zenith", &solzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solzId, NULL, &solzFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solzId, "scale_factor", &solzScale);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solzId, "add_offset", &solzOffset);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "solar_azimuth", &solaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solaId, NULL, &solaFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solaId, "scale_factor", &solaScale);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solaId, "add_offset", &solaOffset);
    check_err(status, __LINE__, __FILE__);

    
    // get IDs for the observations
    status = nc_inq_grp_ncid(ncid_L1B, "observation_data", &observationGrp);
    check_err(status, __LINE__, __FILE__);
    
    // Get varids for each of the Lt_*
    use_rhot = 0;
    status = nc_inq_varid(observationGrp, "Lt_blue", &Lt_blueId);
    if(status == NC_NOERR) {
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(observationGrp, Lt_blueId, NULL, &Lt_blueFillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(observationGrp, "Lt_red", &Lt_redId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(observationGrp, Lt_redId, NULL, &Lt_redFillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(observationGrp, "Lt_SWIR", &Lt_SWIRId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(observationGrp, Lt_SWIRId, NULL, &Lt_SWIRFillValue);
    check_err(status, __LINE__, __FILE__);
    } else { // no Lt, so check for rhot
        status = nc_inq_varid(observationGrp, "rhot_blue", &Lt_blueId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(observationGrp, Lt_blueId, NULL, &Lt_blueFillValue);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_varid(observationGrp, "rhot_red", &Lt_redId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(observationGrp, Lt_redId, NULL, &Lt_redFillValue);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_varid(observationGrp, "rhot_SWIR", &Lt_SWIRId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(observationGrp, Lt_SWIRId, NULL, &Lt_SWIRFillValue);
        check_err(status, __LINE__, __FILE__);
        use_rhot = 1;
    }

    
    file->sd_id = ncid_L1B;
    file->nbands = tot_num_bands;
    file->npix = num_pixels;
    file->nscan = num_scans;
    file->ndets = 1;
    file->terrain_corrected = 1; // presumed.
    file->orbit_number = orbit_number;
    strcpy(file->spatialResolution, "1000 m");

    rdsensorinfo(file->sensorID, l1_input->evalmask,
            "Fobar", (void **) &Fobar);

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
int readl1_ocis(filehandle *file, int32_t line, l1str *l1rec) {

    int status;
    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    //printf("reading oci l1b file\n");
    for (int ip = 0; ip < num_pixels; ip++) {
        l1rec->pixnum[ip] = ip + extract_pixel_start;
    }

    l1rec->npix = file->npix;

    l1rec->scantime = ymds2unix((int16_t)scan_time_year, (int16_t)scan_time_month, (int16_t)scan_time_day, scan_time[line]);

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
    status = nc_get_vara_float(geolocationGrp, latId, start, count, l1rec->lat);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start, count, l1rec->lon);
    check_err(status, __LINE__, __FILE__);

    for(int i=0; i<num_pixels; i++) {
        if(l1rec->lat[i] == latFillValue)
            l1rec->lat[i] = BAD_FLT;
        if(l1rec->lon[i] == lonFillValue)
            l1rec->lon[i] = BAD_FLT;
    }

    status = nc_get_vara_short(geolocationGrp, senaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for(int i=0; i<num_pixels; i++) {
        if(tmpShort[i] == senaFillValue)
            l1rec->sena[i] = BAD_FLT;
        else
            l1rec->sena[i] = tmpShort[i] * senaScale + senaOffset;
    }

    status = nc_get_vara_short(geolocationGrp, senzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for(int i=0; i<num_pixels; i++) {
        if(tmpShort[i] == senzFillValue)
            l1rec->senz[i] = BAD_FLT;
        else
            l1rec->senz[i] = tmpShort[i] * senzScale + senzOffset;
    }

    status = nc_get_vara_short(geolocationGrp, solaId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for(int i=0; i<num_pixels; i++) {
        if(tmpShort[i] == solaFillValue)
            l1rec->sola[i] = BAD_FLT;
        else
            l1rec->sola[i] = tmpShort[i] * solaScale + solaOffset;
    }

    status = nc_get_vara_short(geolocationGrp, solzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for(int i=0; i<num_pixels; i++) {
        if(tmpShort[i] == solzFillValue)
            l1rec->solz[i] = BAD_FLT;
        else
            l1rec->solz[i] = tmpShort[i] * solzScale + solzOffset;
    }
    
    
    start[0] = 0;
    start[1] = line;
    start[2] = 0;
    count[0] = num_blue_bands;
    count[1] = 1;               // 1 line at a time
    count[2] = num_pixels;
    status = nc_get_vara_float(observationGrp, Lt_blueId, start, count, Lt_blue[0]);
    check_err(status, __LINE__, __FILE__);    

    count[0] = num_red_bands;
    status = nc_get_vara_float(observationGrp, Lt_redId, start, count, Lt_red[0]);
    check_err(status, __LINE__, __FILE__);    

    count[0] = num_SWIR_bands;
    status = nc_get_vara_float(observationGrp, Lt_SWIRId, start, count, Lt_SWIR[0]);
    check_err(status, __LINE__, __FILE__);    


    // old
    //Lt = <file val> * <F0 from read_sensor_info> * <esdist()> * cos(solz)

    // this should be the correct equation
    //Lt = <file val> * <F0 from read_sensor_info> * <esdist()>

    
    for(int ip=0; ip<num_pixels; ip++) {
        int band;
        int ib = 0;
        int ipb = ip * tot_num_bands;

        // load up the blue bands skip the last two
        for(band=0; band<expected_num_blue_bands-2; band++) {
            if(Lt_blue[band][ip] == Lt_blueFillValue) {
                l1rec->Lt[ipb] = 0.001; // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1rec->Lt[ipb] = Lt_blue[band][ip];
                if (normalizedLt) {
                    l1rec->Lt[ipb] *= 100;
                } else if (use_rhot) {
                    l1rec->Lt[ipb] *= Fobar[ib] * l1rec->fsol * cos(l1rec->solz[ip]/RADEG) / M_PI;
                } else {
                    l1rec->Lt[ipb] /= 10.; // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the red bands skipping the first two and the last 4
        for(band=2; band<expected_num_red_bands-4; band++) {
            if(Lt_red[band][ip] == Lt_redFillValue) {
                l1rec->Lt[ipb] = 0.001; // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1rec->Lt[ipb] = Lt_red[band][ip];
                if (normalizedLt) {
                    l1rec->Lt[ipb] *= 100;
                } else if (use_rhot) {
                    l1rec->Lt[ipb] *= Fobar[ib] * l1rec->fsol * cos(l1rec->solz[ip]/RADEG) / M_PI;
                } else {
                    l1rec->Lt[ipb] /= 10.; // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the SWIR bands, skip band 3 and 6, hi/low gain wavelengths
        for(band=0; band<expected_num_SWIR_bands; band++) {
            if(band == 3 || band == 6)
                continue;
            if(Lt_SWIR[band][ip] == Lt_SWIRFillValue) {
                l1rec->Lt[ipb] = 0.001; // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1rec->Lt[ipb] = Lt_SWIR[band][ip];
                if (normalizedLt) {
                    l1rec->Lt[ipb] *= 100;
                } else if (use_rhot) {
                    l1rec->Lt[ipb] *= Fobar[ib] * l1rec->fsol * cos(l1rec->solz[ip]/RADEG) / M_PI;
                } else {
                    l1rec->Lt[ipb] /= 10.; // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }
    }
    
    return (LIFE_IS_GOOD);
}


/**
 * Close L1B file, GEO file, and free memory
 * @param file
 * @return 
 */
int closel1_ocis(filehandle *file) {
    int status;

    printf("Closing ocis l1b file\n");
    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);

    // Free memory
    // From openl1_ocis
    if (tmpShort) free(tmpShort);
    if (scan_time) free(scan_time);
    if (Lt_blue) free2d_float(Lt_blue);
    if (Lt_red) free2d_float(Lt_red);
    if (Lt_SWIR) free2d_float(Lt_SWIR);

    return (LIFE_IS_GOOD);
}





