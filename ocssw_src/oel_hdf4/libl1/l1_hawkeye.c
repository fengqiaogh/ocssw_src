/* ============================================================================ */
/* module l1a_hawkeye.c - functions to read HAWKEYE L1A for MSL12               */
/* Written By:  Steve Lockhart SAIC May 2018                                    */
/*      -Started with l1_viirs_nc.c and made changes to support hawkeye.        */
/*                                                                              */
/* ============================================================================ */

/* Issues:
    open:
       -What is spatial resolution? 120m?
    missing fields:
        -Setting orbit number to 0, as it is currently not in L1A.
        -scanQualityId skipped for now
        -Also skipping att angle, so ripples to "Compute polarization rotation angles"
    general:
        -Clean up commented out sections (i.e. the ones that are VIIRS-specific)
        -What are finder_pixels, finder_lines in L1A file?
        -What are dark_pixels in L1A file?
        -Does last line get run twice?
*/

// Includes from l1_viirs_nc.c
#include <nc4utils.h>
#include "l1.h"
#include "libnav.h"
#include <productInfo.h>
#include <allocate2d.h>
#include <allocate3d.h>

// Includes from demo_cal.c
#include <stdio.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <stdlib.h>
#include <fcntl.h>
#include <string.h>
#include <errno.h>
#include <unistd.h>
#include <libgen.h>
#include <math.h>
#include <timeutils.h>
#include <gsl/gsl_math.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort.h>
#include <gsl/gsl_statistics.h>
#include <stdbool.h>

// New includes
#include "l1_hawkeye.h"

static int geoFileId;
static int geoNavigationGrp;
static int geoGeolocationGrp;
static int geoScanLineGrp;
static int l1aScanLineGrp;
static int l1aTelemetryGrp;
static int l1aNavigationGrp;
static int l1aEarthViewGrp;
static int lonId, latId, senzId, senaId, solzId, solaId,
        angId, posId, velId, pixelQualityId;
static int scanDeltaTimeId;
static float *Fobar; // reflectance to radiance conversion factors

static short *tmpShort;
static unsigned char *tmpByte;
static size_t num_scans, num_pixels, num_bands;

static int firstCall = 1;
static double starttime;
static double lastvalidtime;
static int lastvalidscan = 0;
static double time_interval;

static float latGeoFillValue = -999.9;
static float lonGeoFillValue = -999.9;
static short senzGeoFillValue = -32768;
static short senaGeoFillValue = -32768;
static short solzGeoFillValue = -32768;
static short solaGeoFillValue = -32768;

static float latL2FillValue = -999.0;
static float lonL2FillValue = -999.0;
static float senzL2FillValue = -32767;
static float senaL2FillValue = -32767;
static float solzL2FillValue = -32767;
static float solaL2FillValue = -32767;

static int32_t scanDeltaTimeFillValue = -999;
static int32_t scanDeltaTimeValidMin = 0;
static int32_t scanDeltaTimeValidMax = 200000;

// Declarations added for hawkeye L1A.
static size_t num_ccd_temps, num_tlm_blocks, num_tlm_blocks_cleansed;
static double **CCD_temperatures;             // [num_tlm_blocks][num_ccd_temps]
static double **CCD_temperatures_cleansed;    // up to [num_tlm_blocks][num_ccd_temps]
static double CCD_temperature_default = 35.0; // same as K3T[3], used iff all CCD_temperatures are FILL
static double *tlm_delta_time_ms;             // [num_tlm_blocks], ms since start of image
static double *tlm_delta_time_ms_cleansed;    // up to [num_tlm_blocks]
static float fill_value_in_CCD_T;
static int *dn_varid;                         // [num_bands]
static double *CCD_temperatures_this_scan;    // [num_ccd_temps] i.e. interpolated in time
static unsigned short *dn;                   //  [num_pixels], one scan, one band, before cal
static double **bkg_avg;                      // [num_bands][num_pixels], avg dn to subtract
static double **Lt;                           // [num_bands][num_pixels], one scan, after cal
static uint32_t darkSubtracted;
static uint32_t darkHeight;

// radiometric calibration parameters
static size_t cal_num_dates;
static size_t cal_num_temperatures;
static double *time_ref;            // [cal_num_dates]
static double *temperature_ref;     // [cal_num_temperatures]
static double *K1;                  // Gain at nadir pixel, [num_bands]
static double ***K2;                // Temporal correction, BEFORE interpolation, [num_bands][num_pixels][cal_num_dates]
static double **K3;                 // Temperature correction BEFORE interpolation, [num_bands][cal_num_temperatures]
static float **K4;                  // RVS, [num_bands][num_pixels]
static double **K5;                 // Non-linearity (quadratic term), [num_bands][num_pixels]
static float **K7;                  // smile correction [num_bands][num_pixels]
static unsigned short **nominal_dark_counts;

// band co-registration parameters
static int *track_offset;           // [num_bands]
static int *ccd_offset;             // [num_bands]
static float *focal_length;         // [num_bands]
static int reference_band;
static int reference_pixel;

// background subtraction and cropping parameters
static int skip_beglines;
static int skip_endlines;
static int shutter_rampdown;
static int default_darkrows;
static uint32_t firstbad;

// Data cleansing (one-time)
void qc_hawkeye_CCD_T();
// Calibration-related (one-time)
void read_cal_hawkeye(char *cal_path);
// Calibration-related (per scan)
void interp_hawkeye_CCD_T(double delta_time_ms_i);
void calibrate_hawkeye(double current_julian_date, int line);
// Misc utilities
double nan_wmean(double *weights, double *data, size_t n, double fill_value);
void prep_for_interp_double(double *xi, double *x, int N);

/**
 * Open the hawkeye L1A file and perform some one-time tasks (as opposed to tasks that
 * are per scan), including:
 *   -Get L1A dimensions num_scans, num_bands, num_pixels, num_ccd_temps, num_tlm_blocks (static).
 *    Allocate memory for some static arrays based upon these dimensions.
 *   -Get L1A group ids e.g. l1aScanLineGrp, l1aEarthViewGrp and it's 8 "band_%d" var ids (static)
 *   -Get L1A CCD_temperatures and tlm_delta_time_ms, and call qc_hawkeye_CCD_T to create cleansed
 *    versions of these arrays (static).
 *   -Get L1A "time_coverage_start" and "time_coverage_end" to derive time_interval (static) etc.
 *

 * Get
 * @param file
 * @return
 */
int openl1_hawkeye(filehandle *file) {
    char *fltime;

    int ncid_L1A, varid, dimid, status;
    size_t att_len;
    int orbit_number;
    int band_num;
    char nc_search_string[10] = ""; // band_X

    // Open the netcdf4 input file
    if(want_verbose)
        printf("Opening Hawkeye L1A file\n");
    status = nc_open(file->name, NC_NOWRITE, &ncid_L1A);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error opening L1A file \"%s\": %s\n",
                file->name, nc_strerror(status));
        exit(EXIT_FAILURE);
    }
    
    file->private_data = malloc(sizeof(hawkeye_t));
    hawkeye_t *data = (hawkeye_t*)(file->private_data);
    // Get dims from L1A file
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_L1A, "number_of_scans", &dimid));
    nc_inq_dimlen(ncid_L1A, dimid, &num_scans);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_L1A, "number_of_bands", &dimid));
    nc_inq_dimlen(ncid_L1A, dimid, &num_bands);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_L1A, "number_of_pixels", &dimid));
    nc_inq_dimlen(ncid_L1A, dimid, &num_pixels);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_L1A, "ccd_temps", &dimid));
    nc_inq_dimlen(ncid_L1A, dimid, &num_ccd_temps);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_L1A, "number_of_tlm_blocks", &dimid));
    nc_inq_dimlen(ncid_L1A, dimid, &num_tlm_blocks);
    if (want_verbose) {
        printf("Hawkeye L1A npix = %d; nscans = %d; nbands = %d\n",
               (int) num_pixels, (int) num_scans, (int) num_bands);
    }

    // Get hawkeye calibration coefficients from cal file
    // This sets static cal coeffs: K1-K5 as well as time_ref and temperature_ref
    // Do now to get start, end line adjustments
    read_cal_hawkeye(l1_input->calfile);

    // Now that we know dims, prep for additional hawkeye one-time reads by allocating memory
    // for arrays (declared static above). (Memory is freed in closel1_hawkeye.)
    CCD_temperatures = allocate2d_double(num_tlm_blocks, num_ccd_temps);
    CCD_temperatures_cleansed = allocate2d_double(num_tlm_blocks, num_ccd_temps);
    tlm_delta_time_ms = (double *) calloc(num_tlm_blocks, sizeof(double));
    tlm_delta_time_ms_cleansed = (double *) calloc(num_tlm_blocks, sizeof(double));
    dn_varid = (int *) calloc(num_bands, sizeof(int));

    // Get group id from L1A file for GROUP scan_line_attributes.
    if ((nc_inq_grp_ncid(ncid_L1A, "scan_line_attributes", &l1aScanLineGrp))
            != NC_NOERR) {
        fprintf(stderr, "-E- Error finding group \"scan_line_attributes\".\n");
        exit(EXIT_FAILURE);
    }
    TRY_NC(__FILE__, __LINE__,
           nc_inq_varid(l1aScanLineGrp, "delta_time", &scanDeltaTimeId));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(l1aScanLineGrp, scanDeltaTimeId, "_FillValue",
                           &scanDeltaTimeFillValue));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(l1aScanLineGrp, scanDeltaTimeId, "valid_min",
                           &scanDeltaTimeValidMin));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(l1aScanLineGrp, scanDeltaTimeId, "valid_max",
                           &scanDeltaTimeValidMax));

    // Get CCD_temperatures, tlm_delta_time_ms from L1A file from GROUP parameters_telemetry_data
    if ((nc_inq_grp_ncid(ncid_L1A, "parameters_telemetry_data", &l1aTelemetryGrp))
            != NC_NOERR) {
        fprintf(stderr, "-E- Error finding group \"parameters_telemetry_data\".\n");
        exit(EXIT_FAILURE);
    }
    if ((nc_inq_grp_ncid(ncid_L1A, "navigation_data", &l1aNavigationGrp))
            != NC_NOERR) {
        fprintf(stderr, "-E- Error finding group \"navigation_data\" in the L1A file.\n");
        exit(EXIT_FAILURE);
    }
    // Get CCD_temperatures from this GROUP
    TRY_NC(__FILE__, __LINE__,
           nc_inq_varid(l1aTelemetryGrp, "CCD_temperatures", &varid));
    TRY_NC(__FILE__, __LINE__,
           nc_get_var_double(l1aTelemetryGrp, varid, &CCD_temperatures[0][0]));
    if (nc_inq_var_fill(l1aTelemetryGrp, varid, NULL, &fill_value_in_CCD_T)
            != NC_NOERR) {
        fprintf(stderr, "-W- Using default fill value of -999 for CCD_temperatures.\n");
        fill_value_in_CCD_T = -999.0;
    }
    // Get tlm_time_stamp from this GROUP
    TRY_NC(__FILE__, __LINE__,
           nc_inq_varid(l1aTelemetryGrp, "tlm_time_stamp", &varid));
    TRY_NC(__FILE__, __LINE__,
           nc_get_var_double(l1aTelemetryGrp, varid, &tlm_delta_time_ms[0]));

    // was background subtracted on spacecraft?
    if (nc_get_att_uint(l1aTelemetryGrp, NC_GLOBAL, "darkSubtracted", &darkSubtracted)
            != NC_NOERR)
        darkSubtracted = 0;  // assume no, if attribute not present

    // determine number of dark lines at end of image
    if (nc_get_att_uint(l1aTelemetryGrp, NC_GLOBAL, "darkHeight", &darkHeight)
            != NC_NOERR)
        darkHeight = default_darkrows;  // default, if attribute not present

    // Adjust start and end lines
    if (want_verbose) {
        printf("Will skip first %d and last %d scans.\n",
               (int) skip_beglines, (int) darkHeight);
    }
    firstbad = num_scans - darkHeight;

    // Cleanse CCD_temperatures (and corresponding tlm_delta_time_ms), handling fill values
    // This sets static vars CCD_temperatures_cleansed, tlm_delta_time_ms_cleansed, num_tlm_blocks_cleansed
    qc_hawkeye_CCD_T();

    // Get attribute values (from l1_viirs_nc.c) e.g. "time_coverage_start" and "time_coverage_end" to derive
    // time_interval etc. Note that time_coverage_end is the START of the last scan.

    // get start time
    TRY_NC(__FILE__, __LINE__,
           nc_inq_attlen(ncid_L1A, NC_GLOBAL, "time_coverage_start", &att_len));
    fltime = (char *) malloc(att_len + 1); // required space + 1 for trailing null
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_text(ncid_L1A, NC_GLOBAL, "time_coverage_start", fltime));
    fltime[att_len] = '\0';

    // Convert "time_coverage_start" ISO string to unix (seconds since 1/1/1970)
    starttime = lastvalidtime = isodate2unix(fltime);
    free(fltime);

    // get end time
    TRY_NC(__FILE__, __LINE__,
           nc_inq_attlen(ncid_L1A, NC_GLOBAL, "time_coverage_end", &att_len));
    fltime = (char *) malloc(att_len + 1); // required space + 1 for trailing null
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_text(ncid_L1A, NC_GLOBAL, "time_coverage_end", fltime));
    fltime[att_len] = '\0';

    // Convert "time_coverage_stop" ISO string to unix (seconds since 1/1/1970)
    double stoptime = isodate2unix(fltime);
    free(fltime);

    // time_interval may be used in readl1_hawkeye if there is not a good scan time (per scan)
    time_interval = (stoptime - starttime) / (num_scans - 1); // secs per scan

    // orbit number
    if (nc_get_att_int(ncid_L1A, NC_GLOBAL, "orbit_number", &orbit_number) != NC_NOERR)
        orbit_number = 0;
    
    // exposure_ID, roll and time offset
    if (nc_get_att_int(l1aTelemetryGrp, NC_GLOBAL, "exposureID", &data->exposureID) != NC_NOERR)
        data->exposureID = 0;
    if (nc_get_att_float(l1aNavigationGrp, NC_GLOBAL, "time_offset", &data->time_offset) != NC_NOERR)
        data->time_offset = 0.;
    if (nc_get_att_float(l1aNavigationGrp, NC_GLOBAL, "roll_offset", &data->roll_offset) != NC_NOERR)
        data->roll_offset = 0.;

    // Identify the "earth_view_data" GROUP and its "band_%d" vars, to be used later by readl1_hawkeye
    // Store the ids in static variables so we don't have to do nc_inq_grp_ncid per scan.
    if ((nc_inq_grp_ncid(ncid_L1A, "earth_view_data", &l1aEarthViewGrp)) != NC_NOERR) {
        fprintf(stderr, "-E- Error finding group \"earth_view_data\".\n");
        exit(EXIT_FAILURE);
    }
    for (band_num = 0; band_num < num_bands; band_num++) {
        sprintf(nc_search_string, "band_%d", band_num + 1);
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(l1aEarthViewGrp, nc_search_string, &dn_varid[band_num]));
    }

    file->sd_id = ncid_L1A;
    file->nbands = num_bands;
    file->npix = num_pixels;
    file->nscan = num_scans - skip_beglines - darkHeight;
    file->ndets = 1;
    file->terrain_corrected = 0;
    file->orbit_number = orbit_number;
    strcpy(file->spatialResolution, "120 m");

    rdsensorinfo(file->sensorID, l1_input->evalmask,
                 "Fobar", (void **) &Fobar);

    // Setup geofile pointers
    if (file->geofile && file->geofile[0]) {
        if (want_verbose)
            printf("Opening Hawkeye GEO file\n");
        status = nc_open(file->geofile, NC_NOWRITE, &geoFileId);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error opening GEO file \"%s\": %s\n",
                    file->geofile,
                    nc_strerror(status));
            exit(EXIT_FAILURE);
        }

        if ((nc_inq_grp_ncid(geoFileId, "geolocation_data", &geoGeolocationGrp))
                != NC_NOERR) {
            fprintf(stderr, "-E- Error finding group \"geolocation_data\".\n");
            exit(EXIT_FAILURE);
        }
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "longitude", &lonId));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_var_fill(geoGeolocationGrp, lonId, NULL, &lonGeoFillValue));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "latitude", &latId));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_var_fill(geoGeolocationGrp, latId, NULL, &latGeoFillValue));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "sensor_zenith", &senzId));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_var_fill(geoGeolocationGrp, senzId, NULL, &senzGeoFillValue));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "sensor_azimuth", &senaId));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_var_fill(geoGeolocationGrp, senaId, NULL, &senaGeoFillValue));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "solar_zenith", &solzId));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_var_fill(geoGeolocationGrp, solzId, NULL, &solzGeoFillValue));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "solar_azimuth", &solaId));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_var_fill(geoGeolocationGrp, solaId, NULL, &solaGeoFillValue));
        TRY_NC(__FILE__, __LINE__,
               nc_inq_varid(geoGeolocationGrp, "quality_flag", &pixelQualityId));

        if ((nc_inq_grp_ncid(geoFileId, "navigation_data", &geoNavigationGrp))
                != NC_NOERR) {
            fprintf(stderr, "-E- Error finding group \"navigation_data\".\n");
            exit(EXIT_FAILURE);
        }

        TRY_NC(__FILE__, __LINE__, nc_inq_varid(geoNavigationGrp, "att_ang", &angId));
        TRY_NC(__FILE__, __LINE__, nc_inq_varid(geoNavigationGrp, "orb_pos", &posId));
        TRY_NC(__FILE__, __LINE__, nc_inq_varid(geoNavigationGrp, "orb_vel", &velId));
        
        if ((nc_inq_grp_ncid(geoFileId, "scan_line_attributes", &geoScanLineGrp))
                != NC_NOERR) {
            fprintf(stderr, "-E- Error finding group \"scan_line_attributes\".\n");
            exit(EXIT_FAILURE);
        }
        // The following field is not yet in the GEO file
        // TRY_NC(__FILE__, __LINE__,
        //        nc_inq_varid(geoScanLineGrp, "scan_quality", &scanQualityId));
    } // geofile

    // Setup the fill values for the geo products
    productInfo_t* info = allocateProductInfo();
    status = findProductInfo("lat", HAWKEYE, info);
    if (status)
        latL2FillValue = info->fillValue;
    status = findProductInfo("lon", HAWKEYE, info);
    if (status)
        lonL2FillValue = info->fillValue;
    status = findProductInfo("sena", HAWKEYE, info);
    if (status)
        senaL2FillValue = info->fillValue;
    status = findProductInfo("senz", HAWKEYE, info);
    if (status)
        senzL2FillValue = info->fillValue;
    status = findProductInfo("sola", HAWKEYE, info);
    if (status)
        solaL2FillValue = info->fillValue;
    status = findProductInfo("solz", HAWKEYE, info);
    if (status)
        solzL2FillValue = info->fillValue;

    freeProductInfo(info);
    return (LIFE_IS_GOOD);
}

/**
 * Calculate average dark counts for background subtraction
 */
void calc_bkg_hawkeye() {
    int band_num, pixel_num;

    // initialize
    for (pixel_num = 0; pixel_num < num_pixels; pixel_num++)
        for (band_num = 0; band_num < num_bands; band_num++)
            bkg_avg[band_num][pixel_num] = nominal_dark_counts[band_num][pixel_num];

// if dark was not captured, actually earth view, then apply the nominal dark counts in .cal file
    if (!darkSubtracted) {

        // if a full 6000 line scene, calculate dark otherwise assume no valid dark scans
        // and stick with the nominal_dark_counts from the LUT
        if (num_scans == 6000) {
            // find dark lines boundaries
            size_t begline = num_scans - darkHeight + shutter_rampdown;
            size_t endline = num_scans - skip_endlines - 1;
            size_t num_darklines = endline - begline + 1;

            size_t start[] = { 0, 0 };
            size_t count[] = { 1, 1 };
            start[0] = begline;
            count[0] = num_darklines;
            start[1] = 0;
            count[1] = num_pixels;

            // read data for each band
            unsigned short **dn_dark = (unsigned short**) allocate2d_short(num_darklines,
                                                                        num_pixels);
            memset((&dn_dark[0][0]), 0, num_darklines*num_pixels*sizeof(unsigned short));

            for (band_num = 0; band_num < num_bands; band_num++) {
                nc_get_vara_ushort(l1aEarthViewGrp, dn_varid[band_num],
                                start, count, &dn_dark[0][0]);

                // calculate simple average in pixel dimension
                for (pixel_num = 0; pixel_num < num_pixels; pixel_num++) {
                    unsigned long sum_vals = 0;
                    for (size_t line_num = 0; line_num < num_darklines; line_num++)
                        sum_vals += dn_dark[line_num][pixel_num];
                    bkg_avg[band_num][pixel_num] = (double) sum_vals
                            / (double) num_darklines;
                }

            }
            free2d_short((short**) dn_dark);
        }
    }
}

/**
 * Read the specified scan line from the specified L1A file. If this is the first call,
 * call read_cal_hawkeye to read the calibration file, which stores calibration coefficients
 * (for instrument calibration) in static arrays.
 * For each scan, get scan_delta_time_ms and use it to call interp_hawkeye_CCD_T, which sets
 * CCD_temperatures_this_scan[num_ccd_temps]. Then, apply the instrument cal by calling calibrate_hawkeye.
 * Store Lt in l1rec->Lt.
 * Read GEO file.
 *
 * @param file
 * @param line
 * @param l1rec
 * @return
 */
int readl1_hawkeye(filehandle *file, int32_t oline, l1str *l1rec) {

    int i;
    double scan_sec;
    int16_t scan_year, scan_day;
    double Lt_interp;

    size_t start[] = { 0, 0 };
    size_t count[] = { 1, 1 };

    // Additional declarations for hawkeye
    int band_num, pixel_num;
    int16_t scan_month, scan_dom;
    double current_julian_date;
    int32_t scan_delta_time_ms; // ms since start of image per scan

    l1rec->npix = file->npix;
    for (int32_t ip = 0; ip < num_pixels; ip++) {
        l1rec->pixnum[ip] = ip;
    }

    // If first call,
    if (firstCall) {
        firstCall = 0;

        // One-time memory allocations
        tmpShort = (short *) malloc(num_pixels * sizeof(short));
        tmpByte = (unsigned char *) malloc(num_pixels);
        dn = (unsigned short*) calloc(num_pixels, sizeof(unsigned short));
        bkg_avg = allocate2d_double(num_bands, num_pixels);
        Lt = allocate2d_double(num_bands, num_pixels);
        CCD_temperatures_this_scan = (double *) calloc(num_ccd_temps, sizeof(double));

        // Calculate dark counts for background subtraction
        calc_bkg_hawkeye();  // sets bkg_avg

    }
    int32_t line = oline + skip_beglines;

    // Time
    // Get delta_time
    start[0] = line;
    count[0] = 1; // 1 scan at a time
    count[1] = 0;
    TRY_NC(__FILE__, __LINE__,
           nc_get_vara_int(l1aScanLineGrp, scanDeltaTimeId, start, count,
                            &scan_delta_time_ms));
    // Set lastvalidtime, the start of this scan (in secs since 1/1/1970)
    if ((scan_delta_time_ms == scanDeltaTimeFillValue) ||
            (scan_delta_time_ms < scanDeltaTimeValidMin) ||
            (scan_delta_time_ms > scanDeltaTimeValidMax)) {
        l1rec->scantime = lastvalidtime + (time_interval * (line - lastvalidscan));
    } else {
        lastvalidtime = starttime + scan_delta_time_ms / 1000.0;
        lastvalidscan = line;
        l1rec->scantime = lastvalidtime;
    }

    // Set scan_year, scan_day, scan_sec
    unix2yds(l1rec->scantime, &scan_year, &scan_day, &scan_sec);
    // Convert lastvalidtime to julian date
    yd2md(scan_year, scan_day, &scan_month, &scan_dom);
    current_julian_date = jday(scan_year, scan_month, scan_dom)
            + (scan_sec / (24 * 3600.0)) - 0.5; //subtract half a day, as julian time is referenced to noon not midnight

    // Before calibrating this hawkeye scan, we need to interpolate the CCD_temperatures_cleansed (in time).
    // The resulting interpolated array will be CCD_temperatures_this_scan[num_ccd_temps], which is needed by
    // calibrate_hawkeye. Note that the interpolation routine requires the arrays to be of type double.
    interp_hawkeye_CCD_T((double) scan_delta_time_ms);

    /*
       Apply instrument calibration to this scan
       includes:
       - dark count subtraction
       - along-track offset registration
       - application of calibration coefficients
    */
    calibrate_hawkeye(current_julian_date, line);

    for (band_num = 0; band_num < num_bands; band_num++) {

        // select input line according to track offset
        // read only if valid
        int inputline = line + track_offset[band_num];
        if ((-1 < inputline) && (inputline < firstbad)) {
            // apply along-CCD offset and focal length adustments to co-register each pixel
            for (pixel_num = 0; pixel_num < num_pixels; pixel_num++) {

                float inputpix = pixel_num  + ccd_offset[band_num] +
                        ((focal_length[band_num] / focal_length[reference_band] - 1) *
                                (pixel_num - reference_pixel));
                int ipix = floor(inputpix);
                float fpix = inputpix - ipix;

                // if valid, interpolate and store into l1rec->Lt
                if ((-1 < ipix) && (ceil(inputpix) < num_pixels)) {
                    Lt_interp = Lt[band_num][ipix] * (1 - fpix);
                    if (fpix > 0)
                        Lt_interp += Lt[band_num][ipix + 1] * fpix;

                    // ipb = ip * nbands + ib
                    l1rec->Lt[pixel_num * num_bands + band_num] = (float) Lt_interp / 10.0;
                }
            } // pixel_num

        } // valid adjusted input line

    } // band_num

    // If a GEO file was provided, use it...
    if (file->geofile && file->geofile[0]) {

        // set up to read all pixels of the line.
        start[0] = line;
        start[1] = 0;
        count[0] = 1;
        count[1] = num_pixels; // read all pixels

        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_float(geoGeolocationGrp, latId, start, count, l1rec->lat));
        for (i = 0; i < num_pixels; i++)
            if (l1rec->lat[i] == latGeoFillValue)
                l1rec->lat[i] = latL2FillValue;

        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_float(geoGeolocationGrp, lonId, start, count, l1rec->lon));
        for (i = 0; i < num_pixels; i++)
            if (l1rec->lon[i] == lonGeoFillValue)
                l1rec->lon[i] = lonL2FillValue;

        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_short(geoGeolocationGrp, solzId, start, count, tmpShort));
        for (i = 0; i < num_pixels; i++)
            if (tmpShort[i] == solzGeoFillValue)
                l1rec->solz[i] = solzL2FillValue;
            else
                l1rec->solz[i] = tmpShort[i] * 0.01;

        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_short(geoGeolocationGrp, solaId, start, count, tmpShort));
        for (i = 0; i < num_pixels; i++)
            if (tmpShort[i] == solaGeoFillValue)
                l1rec->sola[i] = solaL2FillValue;
            else
                l1rec->sola[i] = tmpShort[i] * 0.01;

        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_short(geoGeolocationGrp, senzId, start, count, tmpShort));
        for (i = 0; i < num_pixels; i++)
            if (tmpShort[i] == senzGeoFillValue)
                l1rec->senz[i] = senzL2FillValue;
            else
                l1rec->senz[i] = tmpShort[i] * 0.01;

        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_short(geoGeolocationGrp, senaId, start, count, tmpShort));
        for (i = 0; i < num_pixels; i++)
            if (tmpShort[i] == senaGeoFillValue)
                l1rec->sena[i] = senaL2FillValue;
            else
                l1rec->sena[i] = tmpShort[i] * 0.01;

        // Load Angles
        float ang[3]; // degrees
        float pos[3]; // km
        float vel[3]; // km/sec
        size_t s3[] = { line, 0 };
        size_t c3[] = { 1, 3 };

        TRY_NC(__FILE__, __LINE__, nc_get_vara_float(geoNavigationGrp, angId, s3, c3, ang));
        TRY_NC(__FILE__, __LINE__, nc_get_vara_float(geoNavigationGrp, posId, s3, c3, pos));
        TRY_NC(__FILE__, __LINE__, nc_get_vara_float(geoNavigationGrp, velId, s3, c3, vel));
        for (i = 0; i < 3; i++) {
            pos[i] /= 1000.; // m   -> km
            vel[i] /= 1000.; // m/s -> km/s
        }

        // Compute polarization rotation angles
        float sen_mat[3][3], coeff[10];
        double mnorm[3];
        ocorient_(pos, vel, ang, sen_mat, coeff);
        for (i = 0; i < 3; i++)
            mnorm[i] = sen_mat[i][0];
        compute_alpha(l1rec->lon, l1rec->lat,
                    l1rec->senz, l1rec->sena,
                    mnorm, l1rec->npix, l1rec->alpha);

        // Check pixel values
        TRY_NC(__FILE__, __LINE__,
            nc_get_vara_uchar(geoGeolocationGrp, pixelQualityId, start, count, tmpByte));
        // 0 Earth_intersection
        // 1 Sensor_zenith
        // 2 Invalid_input
        //
        // make all of the NAVFAIL for now
        for (i = 0; i < num_pixels; i++) {
            if (tmpByte[i])
                l1rec->flags[i] |= NAVFAIL;
        }
    }
    // Earth-sun distance correction for this scan
    int32_t yr = (int32_t) scan_year;
    int32_t dy = (int32_t) scan_day;
    int32_t msec = (int32_t) (scan_sec * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);
    l1rec->fsol = pow(1.0 / esdist, 2);
    l1rec->detnum = 0;

    return (LIFE_IS_GOOD);
}

/**
 * Close L1A file, GEO file, and free memory
 * @param file
 * @return
 */
int closel1_hawkeye(filehandle *file) {

    printf("Closing hawkeye L1A file\n");
    TRY_NC(__FILE__, __LINE__, nc_close(file->sd_id));

    if (file->geofile && file->geofile[0]) {
        printf("Closing Hawkeye GEO file\n");
        TRY_NC(__FILE__, __LINE__, nc_close(geoFileId));

        // Free memory
        // From readl1_hawkeye
        if (tmpShort) free(tmpShort);
        if (tmpByte) free(tmpByte);
        if (CCD_temperatures_this_scan) free(CCD_temperatures_this_scan);
        if (dn) free(dn);
        if (bkg_avg) free2d_double(bkg_avg);
        if (Lt) free2d_double(Lt);
        // From read_cal_hawkeye
        if (temperature_ref) free(temperature_ref);
        if (time_ref) free(time_ref);
        if (K1) free(K1);
        if (K2) free3d_double(K2);
        if (K3) free2d_double(K3);
        if (K4) free2d_float(K4);
        if (K5) free2d_double(K5);
        // From openl1_hawkeye
        if (CCD_temperatures) free2d_double(CCD_temperatures);
        if (CCD_temperatures_cleansed) free2d_double(CCD_temperatures_cleansed);
        if (tlm_delta_time_ms) free(tlm_delta_time_ms);
        if (tlm_delta_time_ms_cleansed) free(tlm_delta_time_ms_cleansed);
        if (dn_varid) free(dn_varid);

    }
    free (file->private_data);
    return (LIFE_IS_GOOD);
}

/**
 * Perform quality control on the CCD_temperatures[num_tlm_blocks][num_ccd_temps]
 * and its time base tlm_delta_time_ms[num_tlm_blocks], creating cleansed versions of both of these arrays.
 * Cleansing includes the following rules:
 *    -If one or more CCD_temperatures for a given tlm block are FILL, replace them with
 *     the mean value of the non-FILL temperatures.
 *    -Identify possibly bad (extreme) values for CCD temperatures e.g. if one
 *     of the 4 CCD thermistors has a value much different from the others.
 *    -If ALL the CCD temperatures for a given tlm block are FILL, skip this tlm block.
 *    -If there are no good tlm blocks, manufacture 2 good rows using a default value for CCD temperature.
 */
void qc_hawkeye_CCD_T() {

    // Misc declarations
    size_t tlm_block_num, tlm_block_cleansed_num, ccd_temp_num;
    bool *tlm_block_good;
    double *CCD_temperatures_Dim2, *nan_mean_weights;
    double CCD_temperature_max_deviation = 10.0;  // deviation from median value
    double *CCD_temperature_sorted;
    size_t *CCD_temperature_sort_index;
    double CCD_temperature_median;

    // Allocate memory to arrays
    tlm_block_good = (bool *) calloc(num_tlm_blocks, sizeof(bool));
    CCD_temperatures_Dim2 = (double *) calloc(num_ccd_temps, sizeof(double));
    nan_mean_weights = (double *) calloc(num_ccd_temps, sizeof(double));
    CCD_temperature_sorted = (double *) calloc(num_ccd_temps, sizeof(double));
    CCD_temperature_sort_index = (size_t *) calloc(num_ccd_temps, sizeof(size_t));

    // If a given row has all FILL values, skip that row.
    for (tlm_block_num = 0; tlm_block_num < num_tlm_blocks; tlm_block_num++) {    // Dim1
        // Apply rules. If at least one of the CCD_temperatures in this tlm_block is good, do not skip this tlm_block.
        tlm_block_good[tlm_block_num] = false;       // Assume worst case
        for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {
            if (CCD_temperatures[tlm_block_num][ccd_temp_num] != fill_value_in_CCD_T) {
                tlm_block_good[tlm_block_num] = true; // At least one good value, so don't skip this row
                break;
            }
        }
    }

    // If one or more CCD_temperatures for a given tlm block are FILL, replace them with the mean value of the
    // non-FILL temperatures.
    // First, set nan_mean_weights
    for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {
        nan_mean_weights[ccd_temp_num] = 1.0;
    }
    tlm_block_cleansed_num = 0;
    for (tlm_block_num = 0; tlm_block_num < num_tlm_blocks; tlm_block_num++) {    // Dim1
        if (tlm_block_good[tlm_block_num]) {
            // This is a good row, so populate CCD_temperatures_cleansed, tlm_delta_time_ms_cleansed and increment
            // tlm_block_cleansed_num
            tlm_delta_time_ms_cleansed[tlm_block_cleansed_num] =
                    tlm_delta_time_ms[tlm_block_num];
            // If any CCD_temperatures are FILL, replace them with the mean value of CCD_temperatures_Dim2
            // We have already removed any rows that have ALL values = FILL, so this should work.
            for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) { // Dim2
                // First pass to build an array just for this row i.e. CCD_temperatures_Dim2
                CCD_temperatures_Dim2[ccd_temp_num] =
                        CCD_temperatures[tlm_block_num][ccd_temp_num];
            }
            for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) { // Dim2
                // Second pass to replace FILL values with the mean for this row
                if (CCD_temperatures_Dim2[ccd_temp_num] == fill_value_in_CCD_T) {
                    CCD_temperatures_cleansed[tlm_block_cleansed_num][ccd_temp_num] =
                            nan_wmean(nan_mean_weights, CCD_temperatures_Dim2,
                                      num_ccd_temps, (double) fill_value_in_CCD_T);
                } else {
                    CCD_temperatures_cleansed[tlm_block_cleansed_num][ccd_temp_num] =
                            CCD_temperatures_Dim2[ccd_temp_num];
                }
            }
            tlm_block_cleansed_num = tlm_block_cleansed_num + 1;
        } else {
            // This is not a good row i.e. all FILL, so do nothing. This row is skipped.
        }
    }

    // At this point, FILL values have been replaced. Identify possibly bad (extreme) values for CCD temperatures e.g. if one
    // of the 4 thermistors has a value much different from the others.
    for (tlm_block_num = 0; tlm_block_num < num_tlm_blocks; tlm_block_num++) {    // Dim1
        for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {    // Dim2
            // First pass to build an array just for this row i.e. CCD_temperatures_Dim2
            CCD_temperatures_Dim2[ccd_temp_num] =
                    CCD_temperatures[tlm_block_num][ccd_temp_num];
        }
        // Sort before calculating median value
        gsl_sort_index(CCD_temperature_sort_index, CCD_temperatures_Dim2, 1,
                       num_ccd_temps);
        for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {
            CCD_temperature_sorted[ccd_temp_num] =
                    CCD_temperatures_Dim2[CCD_temperature_sort_index[ccd_temp_num]];
        }
        // Calculate the median value
        CCD_temperature_median = gsl_stats_median_from_sorted_data(CCD_temperature_sorted,
                                                                   1, num_ccd_temps);
        // Loop through again to see if any values are too far from median value.
        for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {    // Dim2
            if (fabs(CCD_temperatures_Dim2[ccd_temp_num] - CCD_temperature_median)
                    > CCD_temperature_max_deviation) {
                fprintf(stderr,
                        "-W- In tlm block %d, CCD %d records a temperature of %f, more than %f degrees different than the median value of %f\n",
                        (int) tlm_block_num,
                        (int) ccd_temp_num, CCD_temperatures_Dim2[ccd_temp_num],
                        CCD_temperature_max_deviation, CCD_temperature_median);
            }
        }
    }

    // Handle special case where there are no good tlm blocks by setting CCD temperatures to a default value
    if (tlm_block_cleansed_num < 2) {
        // There were NO good tlm blocks
        for (tlm_block_num = 0; tlm_block_num < num_tlm_blocks; tlm_block_num++) { // Dim1
            tlm_delta_time_ms_cleansed[tlm_block_num] = tlm_delta_time_ms[tlm_block_num];
            for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {
                CCD_temperatures_cleansed[tlm_block_num][ccd_temp_num] =
                        CCD_temperature_default;
            }
        }
        num_tlm_blocks_cleansed = num_tlm_blocks;
        fprintf(stderr,
                "-W- There were fewer than 2 good tlm_blocks of CCD_temperatures, so using the default value of %f\n",
                CCD_temperature_default);
    } else {
        // There WAS one or more good tlm blocks
        num_tlm_blocks_cleansed = tlm_block_cleansed_num;
    }

    // Clean up
    free(tlm_block_good);
    free(CCD_temperatures_Dim2);
    free(nan_mean_weights);
    free(CCD_temperature_sorted);
    free(CCD_temperature_sort_index);
}

/**
 * Get calibration coefficients (for instrument calibration) from the specified calibration file
 * and store them in static arrays. These arrays K1-K5, time_ref, and temperature_ref have already been
 * declared static. Here in read_cal_hawkeye, we determine their dimensions, allocate memory to them,
 * and populate them.
 * @param cal_path
 */
void read_cal_hawkeye(char *cal_path) {

    // Declarations for reading cal file
    int status, ncid_CAL, varid, dimid;
    int calGrp, geoGrp;

    // Misc declarations
    size_t cal_num_bands = 0;
    size_t cal_num_pixels = 0;

    // Open LUT for calibrations
    printf("Reading Hawkeye calibration LUT\n");
    status = nc_open(cal_path, NC_NOWRITE, &ncid_CAL);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error opening CAL file \"%s\": %s\n",
                cal_path, nc_strerror(status));
        exit(EXIT_FAILURE);
    }

    // Get dims from cal file
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_CAL, "number_of_bands", &dimid));
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_bands);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_CAL, "number_of_pixels", &dimid));
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_pixels);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_CAL, "number_of_times", &dimid));
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_dates);
    TRY_NC(__FILE__, __LINE__, nc_inq_dimid(ncid_CAL, "number_of_temperatures", &dimid));
    nc_inq_dimlen(ncid_CAL, dimid, &cal_num_temperatures);

    // Check to see if cal file dims = L1A file dims
    if (cal_num_bands != num_bands) {
        fprintf(stderr, "-E- num_bands in cal file not equal to num_bands in L1A file\n");
        exit(EXIT_FAILURE);
    }
    if (cal_num_pixels != num_pixels) {
        fprintf(stderr,
                "-E- num_pixels in cal file not equal to num_pixels in L1A file\n");
        exit(EXIT_FAILURE);
    }

    // Allocate space for arrays used for calibration coefficients
    // 1D arrays
    time_ref = (double *) calloc(cal_num_dates, sizeof(double));
    temperature_ref = (double *) calloc(cal_num_temperatures, sizeof(double));
    K1 = (double *) calloc(num_bands, sizeof(double));
    // 2D arrays
    K3 = allocate2d_double(num_bands, cal_num_temperatures);
    K4 = allocate2d_float(num_bands, num_pixels);
    K5 = allocate2d_double(num_bands, num_pixels);
    // 3D arrays
    K2 = allocate3d_double(num_bands, num_pixels, cal_num_dates);
    nominal_dark_counts = (unsigned short**)  allocate2d_short(num_bands, num_pixels);

    // Get variables from group radiometric_calibration
    if ((nc_inq_grp_ncid(ncid_CAL, "radiometric_calibration", &calGrp)) != NC_NOERR) {
        fprintf(stderr, "-E- Error finding group \"radiometric_calibration\".\n");
        exit(EXIT_FAILURE);
    }
    // Get K1
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K1", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_double(calGrp, varid, K1));
    // Get K2
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K2", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_double(calGrp, varid, &K2[0][0][0]));
    // Get K2t
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K2t", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_double(calGrp, varid, time_ref));
    // Get K3
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K3", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_double(calGrp, varid, &K3[0][0]));
    // Get K3T
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K3T", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_double(calGrp, varid, temperature_ref));
    // Get K4
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K4", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_float(calGrp, varid, &K4[0][0]));
    // Get K5
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K5", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_double(calGrp, varid, &K5[0][0]));

    // Get K7
    if (l1_input->rad_opt) {
        TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "K7", &varid));
        K7 = allocate2d_float(num_bands, num_pixels);
        TRY_NC(__FILE__, __LINE__, nc_get_var_float(calGrp, varid, &K7[0][0]));
    }
    // Get nominal_dark_counts
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(calGrp, "nominal_dark_counts", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_ushort(calGrp, varid, &nominal_dark_counts[0][0]));

    // Get variables from group geometric_parameters
    if ((nc_inq_grp_ncid(ncid_CAL, "geometric_parameters", &geoGrp)) != NC_NOERR) {
        fprintf(stderr, "-E- Error finding group \"geometric_parameters\".\n");
        exit(EXIT_FAILURE);
    }
    // band co-registration parameters
    track_offset = calloc(num_bands, sizeof(int));
    ccd_offset = calloc(num_bands, sizeof(int));
    focal_length = (float *) calloc(num_bands, sizeof(float));
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(geoGrp, "track_offset", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_int(geoGrp, varid, track_offset));
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(geoGrp, "CCD_offset", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_int(geoGrp, varid, ccd_offset));
    TRY_NC(__FILE__, __LINE__, nc_inq_varid(geoGrp, "focal_length", &varid));
    TRY_NC(__FILE__, __LINE__, nc_get_var_float(geoGrp, varid, focal_length));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(geoGrp, NC_GLOBAL, "reference_band", &reference_band));
    reference_band -= 1; // set to zero-based
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(geoGrp, NC_GLOBAL, "reference_pixel", &reference_pixel));
    reference_pixel -= 1; // set to zero-based 
    // background subtraction and cropping parameters
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(geoGrp, NC_GLOBAL, "skip_beglines", &skip_beglines));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(geoGrp, NC_GLOBAL, "skip_endlines", &skip_endlines));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(geoGrp, NC_GLOBAL, "shutter_rampdown", &shutter_rampdown));
    TRY_NC(__FILE__, __LINE__,
           nc_get_att_int(geoGrp, NC_GLOBAL, "default_darkrows", &default_darkrows));

    // Close file
    TRY_NC(__FILE__, __LINE__, nc_close(ncid_CAL));

}

/**
 * The array CCD_temperatures_cleansed has dimensions num_tlm_blocks_cleansed x num_ccd_temps.
 * Given the input value of delta_time_ms_i, we interpolate on the first dimension, generating
 * the array CCD_temperatures_this_scan, having dimension 1 x num_ccd_temps
 * @param delta_time_ms_i
 */
void interp_hawkeye_CCD_T(double delta_time_ms_i) {

    // Misc declarations
    size_t tlm_block_num, ccd_temp_num;
    double *CCD_temperatures_Dim1;

    // Allocate memory for arrays
    CCD_temperatures_Dim1 = (double *) calloc(num_tlm_blocks_cleansed, sizeof(double));

    gsl_interp_accel *acc_delta_time = gsl_interp_accel_alloc();
    gsl_spline *spline_delta_time = gsl_spline_alloc(gsl_interp_linear,
                                                     num_tlm_blocks_cleansed);

    // First, before interpolating, make sure delta_time_ms_i is not outside the tlm_delta_time_ms_cleansed array
    prep_for_interp_double(&delta_time_ms_i, tlm_delta_time_ms_cleansed,
                           num_tlm_blocks_cleansed);
    for (ccd_temp_num = 0; ccd_temp_num < num_ccd_temps; ccd_temp_num++) {      // Dim2
        for (tlm_block_num = 0; tlm_block_num < num_tlm_blocks_cleansed;
                tlm_block_num++) {                                              // Dim1
            // Temporarily store a column in CCD_temperatures_Dim1
            CCD_temperatures_Dim1[tlm_block_num] =
                    CCD_temperatures_cleansed[tlm_block_num][ccd_temp_num];
        }
        // Interpolate this column in time to get CCD_temperatures_i[ccd_temp_num]
        gsl_spline_init(spline_delta_time, tlm_delta_time_ms_cleansed,
                        &CCD_temperatures_Dim1[0], num_tlm_blocks_cleansed);
        CCD_temperatures_this_scan[ccd_temp_num] = gsl_spline_eval(spline_delta_time,
                                                                   delta_time_ms_i,
                                                                   acc_delta_time);
    }

    // Clean up
    free(CCD_temperatures_Dim1);
    gsl_spline_free(spline_delta_time);
    gsl_interp_accel_free(acc_delta_time);
}

/**
 * Apply instrument calibration, converting dn into Lt.
 * @param current_julian_date
 * @return
 */
void calibrate_hawkeye(double current_julian_date, int line) {

    double this_dn;
    double current_temp_C;

    size_t start[] = { 0, 0 };
    size_t count[] = { 1, 1 };

    int32_t band_num, pixel_num, ccd_temp_num;
    double cal_coeff;

    // Declarations for interpolation
    double **K2i; // Temporal correction AFTER interpolation, [num_bands][num_pixels]
    double *K3i;  // Temperature correction AFTER interpolation, [num_bands]
    gsl_interp_accel *acc_T = gsl_interp_accel_alloc();
    gsl_interp_accel *acc_date = gsl_interp_accel_alloc();
    gsl_spline *spline_T = gsl_spline_alloc(gsl_interp_linear, cal_num_temperatures);
    //gsl_spline *spline_date = gsl_spline_alloc(gsl_interp_linear, cal_num_dates);

    // Allocate memory for arrays
    K2i = allocate2d_double(num_bands, num_pixels);
    K3i = (double *) calloc(num_bands, sizeof(double));

    // The original intent of K2 was to account for a temporal decay.
    // This is being overridden here to apply a simple time based
    // incremental (stepwise) adjustment to K4.  Should a decay be necessary,
    // this will need to be rethought...
    int32_t K2t_idx;
    for (K2t_idx = 0; K2t_idx < cal_num_dates; K2t_idx++){
        // find the index where current time is greater than the index time
        if ((time_ref[K2t_idx] - current_julian_date) >= 0)
            break;
    }
    K2t_idx--;

    // Interpolate in time for K2, picking some time just for demo.
    prep_for_interp_double(&current_julian_date, time_ref, cal_num_dates);

    for (band_num = 0; band_num < num_bands; band_num++) {
        for (pixel_num = 0; pixel_num < num_pixels; pixel_num++) {
            K2i[band_num][pixel_num] = K2[band_num][pixel_num][K2t_idx];
            // See note above about change to K2 behavior  This loop left in just in case
            // (this could have been accomplished without K2i, just index K2 by K2t_idx
            // in the calibration loop below...)

            // gsl_spline_init(spline_date, time_ref, &K2[band_num][pixel_num][0],
            //                 cal_num_dates);
            // K2i[band_num][pixel_num] = gsl_spline_eval(spline_date, current_julian_date,
            //                                           acc_date);
        }
    }

    // Interpolate in temperature for K3
    for (band_num = 0; band_num < num_bands; band_num++) {
        // Given the band_num, pick the correct CCD_temperature
        ccd_temp_num = (int) (band_num / 2);
        current_temp_C = CCD_temperatures_this_scan[ccd_temp_num];
        // Make sure this temperature is in bounds
        prep_for_interp_double(&current_temp_C, temperature_ref, cal_num_temperatures);
        // Interpolate K3, given this temperature
        gsl_spline_init(spline_T, temperature_ref, &K3[band_num][0],
                        cal_num_temperatures);
        K3i[band_num] = gsl_spline_eval(spline_T, current_temp_C, acc_T);
    }

    // Apply cal to sample scan
    count[0] = 1; // 1 line at a time
    count[1] = num_pixels;
    for (band_num = 0; band_num < num_bands; band_num++) {
        // select input line according to track offset
        // read only if valid
        int inputline = line + track_offset[band_num];
        if ((-1 < inputline) && (inputline < firstbad)) {
            start[0] = inputline;
            TRY_NC(__FILE__, __LINE__,
                nc_get_vara_ushort(l1aEarthViewGrp, dn_varid[band_num],
                                    start, count, &dn[0]));

            for (pixel_num = 0; pixel_num < num_pixels; pixel_num++) {
                // Dark subtraction
                this_dn = (double) dn[pixel_num] - bkg_avg[band_num][pixel_num];
                /* if (this_dn < (-1 * bkg_avg[band_num][pixel_num])) */
                /*     printf("B%d\tP%d\tdn - bkg_avg = %f\n", band_num, pixel_num, this_dn); */
                if (this_dn < 0)
                    this_dn = 0;

                cal_coeff = (K1[band_num]) * (K2i[band_num][pixel_num]) * (K3i[band_num])
                        * (K4[band_num][pixel_num]);

                // apply smile correction
                if (l1_input->rad_opt)
                    cal_coeff *= K7[band_num][pixel_num];

                Lt[band_num][pixel_num] = this_dn * cal_coeff
                        * (1 + (K5[band_num][pixel_num]) * this_dn);
            }
        }
    }

    // Clean up e.g. free memory
    free2d_double(K2i);
    free(K3i);
    gsl_spline_free(spline_T);
    gsl_interp_accel_free(acc_T);
    //gsl_spline_free(spline_date);
    gsl_interp_accel_free(acc_date);

}

double nan_wmean(double *weights, double *data, size_t n, double fill_value) {
    size_t i, fill_count;
    double wmean;
    // Set weight to zero if corresponding data element is FILL
    fill_count = 0;
    for (i = 0; i < n; i++) {
        if (data[i] == fill_value) {
            weights[i] = 0;
            fill_count++;
        }
    }
    // Call gsl routine for weighted mean
    if (fill_count < n) {
        wmean = gsl_stats_wmean(weights, 1, data, 1, n);
    } else {
        wmean = fill_value;
    }
    return wmean;

}

/**
 * If xi is outside the range of x, set it to the nearest neighbor.
 * @param xi
 * @param x
 * @param N
 */
void prep_for_interp_double(double *xi, double *x, int N) {
    // Adjust xi to be within range of x
    double xmin, xmax;
    gsl_stats_minmax(&xmin, &xmax, x, 1, N);
    if (*xi < xmin) {
        *xi = xmin;
    }
    if (*xi > xmax) {
        *xi = xmax;
    }
}
