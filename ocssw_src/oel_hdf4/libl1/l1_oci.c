/* ============================================================================ */
/* module l1b_oci.c - functions to read OCI L1B for l2gen           */
/* Written By:  Don Shea                                                        */
/*                                                                              */
/* ============================================================================ */

/* Issues:
        -Assume num_pixels = ccd_pixels = SWIR_pixels
*/

/*
 further var population by Tincho, 6/16->
*/ 

#include <netcdf.h>
#include "l1_oci.h"

#include "l1.h"
#include <nc4utils.h>
#include "libnav.h"
#include <stdio.h>
#include <math.h>
#include <allocate2d.h>
#include "l1_oci_private.h"

#define SATURATION_BIT 1

//static float *Fobar; // reflectance to radiance conversion factors
static int extract_pixel_start = 0;
static int use_rhot = 0;

static int bad_num_bands = 0;

static short *tmpShort;

// whole file stuff
static size_t expected_num_blue_bands = 119;
static size_t expected_num_red_bands = 163;
static size_t expected_num_SWIR_bands = 9;

static size_t num_scans, num_pixels;
static size_t num_blue_bands, num_red_bands, num_SWIR_bands;
static size_t tot_num_bands = 286;
static int ncid_L1B;

// scan line attributes
static double *scan_time; // seconds of day
static double file_start_day; // unix time of start day 00:00:00
static float *tilt;
static unsigned char *scanQual;
static uint8_t *hamside;

// geolocation data
static int geolocationGrp; // netCDF groupid
static int lonId, latId, heightId, senzId, senaId, solzId, solaId; // netCDF varids
static float latFillValue = BAD_FLT;
static float lonFillValue = BAD_FLT;
static short heightFillValue = BAD_FLT;
static short senzFillValue = BAD_FLT;
static float senzScale = 0.01;
static float senzOffset = 0.0;
static short senaFillValue = BAD_FLT;
static float senaScale = 0.01;
static float senaOffset = 0.0;
static short solzFillValue = BAD_FLT;
static float solzScale = 0.01;
static float solzOffset = 0.0;
static short solaFillValue = BAD_FLT;
static float solaScale = 0.01;
static float solaOffset = 0.0;

// Observation data
static int observationGrp,navigationGrp;
static int Lt_blueId, Lt_redId, Lt_SWIRId, tiltId, qual_SWIRId;
static int ccdScanAnglesId;
static float **Lt_blue;                //[num_blue_bands][num_pixels], This scan
static float **Lt_red;                 //[num_red_bands][num_pixels], This scan
static float **Lt_SWIR;                //[num_SWIR_bands][num_pixels], This scan
static float Lt_blueFillValue = BAD_FLT;
static float Lt_redFillValue = BAD_FLT;
static float Lt_SWIRFillValue = BAD_FLT;
static float tiltFillValue = BAD_FLT;
static uint8_t **qual_SWIR;
static float *blue_solar_irradiance; // [num_blue_bands]
static float *red_solar_irradiance; // [num_red_bands]
static float *SWIR_solar_irradiance; // [num_SWIR_bands]
static double earth_sun_distance_correction;

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
int openl1_oci(filehandle * file) {
    int dimid, status;
    
    // Open the netcdf4 input file
    status = nc_open(file->name, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n",
                __FILE__, __LINE__, file->name);
        return (1);
    }
    
    // num_scans
    status = nc_inq_dimid(ncid_L1B, "scans", &dimid);
    if (status != NC_NOERR) {
        // Try to look for the old name for this dimension
        if (nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid) != NC_NOERR) {
            fprintf(stderr, "-E- Error reading scan dimension.\n");
            exit(EXIT_FAILURE);
        }
    }
    nc_inq_dimlen(ncid_L1B, dimid, &num_scans);

    // num_pixels
    status = nc_inq_dimid(ncid_L1B, "pixels", &dimid);
    if (status != NC_NOERR) {
        // If "pixels" doesn't exist, it's the old "ccd_pixels"
        if (nc_inq_dimid(ncid_L1B, "ccd_pixels", &dimid)) {
            fprintf(stderr, "-E- Error reading num_pixels.\n");
            exit(EXIT_FAILURE);
        }
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
        fprintf(stderr, "-W- Not enough blue bands, expecting %d, found %d.\n",
                (int)expected_num_blue_bands, (int)num_blue_bands);
        bad_num_bands = 1;
    }

    // num_red_bands
    status = nc_inq_dimid(ncid_L1B, "red_bands", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_red_bands.\n");
        exit(EXIT_FAILURE);
    };
    nc_inq_dimlen(ncid_L1B, dimid, &num_red_bands);
    if(num_red_bands < expected_num_red_bands) {
        fprintf(stderr, "-W- Not enough red bands, expecting %d, found %d.\n",
                (int)expected_num_red_bands, (int)num_red_bands);
        bad_num_bands = 1;
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
        printf("OCI L1B Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);
    } // want_verbose

    // allocate all of the data
    tmpShort = (short*) malloc(num_pixels * sizeof(short));
    scan_time = (double*) malloc(num_scans * sizeof(double));

    Lt_blue = allocate2d_float(num_blue_bands, num_pixels);
    Lt_red = allocate2d_float(num_red_bands, num_pixels);
    Lt_SWIR = allocate2d_float(num_SWIR_bands, num_pixels);
    tilt = (float*) malloc(num_scans * sizeof(float));
    hamside = (uint8_t *) malloc(num_scans * sizeof(uint8_t));
    scanQual = (unsigned char *) malloc(num_scans * sizeof(unsigned char));
    qual_SWIR = allocate2d_uchar(num_SWIR_bands, num_pixels);

    // Get group id from L1B file for GROUP scan_line_attributes.
    int scanLineGrp;
    if ((nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &scanLineGrp)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding scan_line_attributes.\n");
        exit(EXIT_FAILURE);
    }   
    int varId;
    double scan_timeFillValue = BAD_FLT;
    status = nc_inq_varid(scanLineGrp, "time", &varId);
    if(status == NC_NOERR) {
        status = nc_inq_var_fill(scanLineGrp, varId, NULL, &scan_timeFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(scanLineGrp, varId, scan_time);
        check_err(status, __LINE__, __FILE__);    
    } else {
        status = nc_inq_varid(scanLineGrp, "ev_mid_time", &varId);
        check_err(status, __LINE__, __FILE__);    
        status = nc_inq_var_fill(scanLineGrp, varId, NULL, &scan_timeFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(scanLineGrp, varId, scan_time);
        check_err(status, __LINE__, __FILE__);    
    }

    // get HAM side info
    status = nc_inq_varid(scanLineGrp, "HAM_side", &varId);
    if( ( status = nc_inq_varid(scanLineGrp, "HAM_side", &varId) ) == NC_NOERR) {
        status = nc_get_var_ubyte(scanLineGrp, varId, hamside);
        check_err(status, __LINE__, __FILE__);
    } else {
        fprintf(stderr, "-E- Error reading the HAM side data.\n");
        exit(EXIT_FAILURE);
    }

    // get scan quality info
    status = nc_inq_varid(scanLineGrp, "scan_quality_flags", &varId);
    if( ( status = nc_inq_varid(scanLineGrp, "scan_quality_flags", &varId) ) == NC_NOERR) {
        status = nc_get_var_uchar(scanLineGrp, varId, scanQual);
        check_err(status, __LINE__, __FILE__);
    } else {
        fprintf(stderr, "-E- Error reading the scan quality flags.\n");
        exit(EXIT_FAILURE);
    }

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
    file_start_day = ymds2unix(syear, smon, sday, 0.0);

    free(time_str);

    for(int i=0; i<num_scans; i++) {
        if(scan_time[i] == scan_timeFillValue)
            scan_time[i] = BAD_FLT;
    }
    
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
    status = nc_inq_varid(geolocationGrp, "height", &heightId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, heightId, NULL, &heightFillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "sensor_zenith", &senzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, senzId, NULL, &senzFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senzId, "scale_factor", &senzScale);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senzId, "add_offset", &senzOffset);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);
    
    status = nc_inq_varid(geolocationGrp, "sensor_azimuth", &senaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, senaId, NULL, &senaFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senaId, "scale_factor", &senaScale);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, senaId, "add_offset", &senaOffset);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "solar_zenith", &solzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solzId, NULL, &solzFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solzId, "scale_factor", &solzScale);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solzId, "add_offset", &solzOffset);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "solar_azimuth", &solaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solaId, NULL, &solaFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solaId, "scale_factor", &solaScale);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solaId, "add_offset", &solaOffset);
    // ignore missing
    //check_err(status, __LINE__, __FILE__);

    
    // get IDs for the observations
    status = nc_inq_grp_ncid(ncid_L1B, "observation_data", &observationGrp);
    check_err(status, __LINE__, __FILE__);
    
    // Get varids for each of the Lt_*
    status = nc_inq_varid(observationGrp, "Lt_blue", &Lt_blueId);
    if(status == NC_NOERR) {
        use_rhot = 0;

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
        use_rhot = 1;

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
    }

    status = nc_inq_grp_ncid(ncid_L1B, "navigation_data", &navigationGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(navigationGrp, "tilt_angle", &tiltId);
    if(status){
        status = nc_inq_varid(navigationGrp, "tilt", &tiltId);
        if(status)
            check_err(status, __LINE__, __FILE__);
    }
    status = nc_inq_var_fill(navigationGrp, tiltId, NULL, &tiltFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(navigationGrp,tiltId,tilt);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(navigationGrp, "CCD_scan_angles", &ccdScanAnglesId);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(observationGrp, "qual_SWIR", &qual_SWIRId);
    check_err(status, __LINE__, __FILE__);

    if(use_rhot) {
        //rdsensorinfo(file->sensorID, l1_input->evalmask, "Fobar", (void **) &Fobar);

        status = nc_get_att_double(ncid_L1B, NC_GLOBAL, "earth_sun_distance_correction", &earth_sun_distance_correction);
        check_err(status, __LINE__, __FILE__);

        int sensorBandGrp;
        int tmpId;
        blue_solar_irradiance = malloc(num_blue_bands * sizeof(float));
        red_solar_irradiance = malloc(num_red_bands * sizeof(float));
        SWIR_solar_irradiance = malloc(num_SWIR_bands * sizeof(float));
        status = nc_inq_grp_ncid(ncid_L1B, "sensor_band_parameters", &sensorBandGrp);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_varid(sensorBandGrp, "blue_solar_irradiance", &tmpId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(sensorBandGrp, tmpId, blue_solar_irradiance);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_varid(sensorBandGrp, "red_solar_irradiance", &tmpId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(sensorBandGrp, tmpId, red_solar_irradiance);
        check_err(status, __LINE__, __FILE__);

        status = nc_inq_varid(sensorBandGrp, "SWIR_solar_irradiance", &tmpId);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(sensorBandGrp, tmpId, SWIR_solar_irradiance);
        check_err(status, __LINE__, __FILE__);
    }

    file->sd_id = ncid_L1B;
    file->nbands = tot_num_bands;
    file->npix = num_pixels;
    file->nscan = num_scans;
    file->ndets = 1;
    file->terrain_corrected = 1; // presumed.
    strcpy(file->spatialResolution, "1000 m");

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
 * @param lonlat if 1, reads only the data needed for lon and lat.
 * @return 
 * 
 */
int readl1_oci(filehandle *file, int32_t line, l1str *l1rec, int lonlat) {
    int status;

    PolcorOciData* polcorPrivateData = (PolcorOciData*)l1rec->private_data;

    // allocate private data for OCI that is needed for polarization
    // dont run if running lonlat mode
    if (l1rec->private_data == NULL && !lonlat) {
        polcorPrivateData = (PolcorOciData*)malloc(sizeof (PolcorOciData));
        polcorPrivateData->ccdScanAngle = (float *)malloc(num_pixels * sizeof(float));

        // read bands for m12 and m13 coefs
        polcorPrivateData->blueM12Coefs = (float *)malloc((num_blue_bands * HAM_SIDES * POLARIZATION_COEFS) * sizeof(float));
        polcorPrivateData->blueM13Coefs = (float *)malloc((num_blue_bands * HAM_SIDES * POLARIZATION_COEFS) * sizeof(float));

        // red band for m12 and m13 coefs
        polcorPrivateData->redM12Coefs = (float *)malloc((num_red_bands * HAM_SIDES * POLARIZATION_COEFS) * sizeof(float));
        polcorPrivateData->redM13Coefs = (float *)malloc((num_red_bands * HAM_SIDES * POLARIZATION_COEFS) * sizeof(float));

        // swir band for m12 and m13 coefs
        polcorPrivateData->swirM12Coefs = (float *)malloc((num_SWIR_bands * HAM_SIDES * POLARIZATION_COEFS) * sizeof(float));
        polcorPrivateData->swirM13Coefs = (float *)malloc((num_SWIR_bands * HAM_SIDES * POLARIZATION_COEFS) * sizeof(float));

        // hook it onto l1rec's private_data
        l1rec->private_data = polcorPrivateData;

        // save the bands on initial read because it wont change for the rest of the scan
        polcorPrivateData->num_blue_bands = num_blue_bands;
        polcorPrivateData->num_red_bands = num_red_bands;
        polcorPrivateData->num_swir_bands = num_SWIR_bands;

        // read static data

        // grab sensor_band_parameters group and ids for bands M12 an M13 coefs
        int sensorBandGroup = -1;
        int blueM12Id = -1;
        int blueM13Id = -1;
        int redM12Id = -1;
        int redM13Id = -1;
        int swirM12Id = -1;
        int swirM13Id = -1;

        status = nc_inq_grp_ncid(ncid_L1B, "sensor_band_parameters", &sensorBandGroup);
        check_err(status, __LINE__, __FILE__);

        // grab the ids for band m12 and m13
        status = nc_inq_varid(sensorBandGroup, "blue_m12_coef", &blueM12Id);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(sensorBandGroup, "blue_m13_coef", &blueM13Id);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(sensorBandGroup, "red_m12_coef", &redM12Id);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(sensorBandGroup, "red_m13_coef", &redM13Id);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(sensorBandGroup, "SWIR_m12_coef", &swirM12Id);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_varid(sensorBandGroup, "SWIR_m13_coef", &swirM13Id);
        check_err(status, __LINE__, __FILE__);

        // grabbing blue M12 and M13
        status = nc_get_var_float(sensorBandGroup, blueM12Id, polcorPrivateData->blueM12Coefs);
        check_err(status, __LINE__, __FILE__);
        status =nc_get_var_float(sensorBandGroup, blueM13Id, polcorPrivateData->blueM13Coefs);
        check_err(status, __LINE__, __FILE__);

        // grabbing red M12 and M13
        status =nc_get_var_float(sensorBandGroup, redM12Id, polcorPrivateData->redM12Coefs);
        check_err(status, __LINE__, __FILE__);
        status =nc_get_var_float(sensorBandGroup, redM13Id, polcorPrivateData->redM13Coefs);
        check_err(status, __LINE__, __FILE__);
        
        // grabbing swir M12 and M13
        status =nc_get_var_float(sensorBandGroup, swirM12Id, polcorPrivateData->swirM12Coefs);
        check_err(status, __LINE__, __FILE__);
        status =nc_get_var_float(sensorBandGroup, swirM13Id, polcorPrivateData->swirM13Coefs);
        check_err(status, __LINE__, __FILE__);

    }

    size_t start[] = { 0, 0, 0 };
    size_t count[] = { 1, 1, 1 };

    for (int ip = 0; ip < num_pixels; ip++) {
        l1rec->pixnum[ip] = ip + extract_pixel_start;
    }

    l1rec->npix = file->npix;
    l1rec->scantime = BAD_FLT;
    l1rec->mside = hamside[line];
    if (scan_time[line] != BAD_FLT) {
        l1rec->scantime = file_start_day + scan_time[line];
    }

    l1rec->tilt=tilt[line];

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

    // dont need height for lonlat, but do when it is not lonlat
    if (!lonlat) {
        status = nc_get_vara_float(geolocationGrp, heightId, start, count, l1rec->height);
        check_err(status, __LINE__, __FILE__);
    }

    for(int i=0; i<num_pixels; i++) {
        if(l1rec->lat[i] == latFillValue)
            l1rec->lat[i] = BAD_FLT;
        if(l1rec->lon[i] == lonFillValue)
            l1rec->lon[i] = BAD_FLT;
    }

    if(scanQual[line] & 1) {
        for(int i=0; i<num_pixels; i++)
            l1rec->navwarn[i] = 1;
    }

    status = nc_get_vara_short(geolocationGrp, solzId, start, count, tmpShort);
    check_err(status, __LINE__, __FILE__);
    for(int i=0; i<num_pixels; i++) {
        if(tmpShort[i] == solzFillValue)
            l1rec->solz[i] = BAD_FLT;
        else
            l1rec->solz[i] = tmpShort[i] * solzScale + solzOffset;
    }

    // lon lat, only need solarz and nothing else
    // at this point, it should have read in: time, lon, lat, solz
    if (lonlat)
        return (LIFE_IS_GOOD);

    if (bad_num_bands) {
        fprintf(stderr, "-E- Bad number of bands\n");
        exit(EXIT_FAILURE);
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
    
    
    start[0] = 0;
    start[1] = line;
    start[2] = 0;
    count[0] = expected_num_blue_bands;
    count[1] = 1;               // 1 line at a time
    count[2] = num_pixels;
    status = nc_get_vara_float(observationGrp, Lt_blueId, start, count, Lt_blue[0]);
    check_err(status, __LINE__, __FILE__);    

    count[0] = expected_num_red_bands;
    status = nc_get_vara_float(observationGrp, Lt_redId, start, count, Lt_red[0]);
    check_err(status, __LINE__, __FILE__);    

    count[0] = expected_num_SWIR_bands;
    status = nc_get_vara_float(observationGrp, Lt_SWIRId, start, count, Lt_SWIR[0]);
    check_err(status, __LINE__, __FILE__);    

    status = nc_get_vara_ubyte(observationGrp, qual_SWIRId, start, count, qual_SWIR[0]);
    check_err(status, __LINE__, __FILE__);
    
    for(int ip=0; ip<num_pixels; ip++) {
        int band;
        int ib = 0;
        int ipb = ip * tot_num_bands;

        // load up the blue bands skip last 3
        for(band=0; band<expected_num_blue_bands-3; band++) {
            if(Lt_blue[band][ip] == Lt_blueFillValue) {
                l1rec->Lt[ipb] = BAD_FLT;
            } else {
                l1rec->Lt[ipb] = Lt_blue[band][ip];
                if (use_rhot) {
                    //l1rec->Lt[ipb] *= Fobar[ib] * l1rec->fsol * cos(l1rec->solz[ip]/OEL_RADEG) / OEL_PI;
                    l1rec->Lt[ipb] *= blue_solar_irradiance[band] * cos(l1rec->solz[ip]/OEL_RADEG) / earth_sun_distance_correction / OEL_PI / 10.0;
                } else {
                    l1rec->Lt[ipb] /= 10.; // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the red bands use all of them
        for(band=0; band<expected_num_red_bands; band++) {
            if(Lt_red[band][ip] == Lt_redFillValue) {
                l1rec->Lt[ipb] = BAD_FLT;
            } else {
                l1rec->Lt[ipb] = Lt_red[band][ip];
                if (use_rhot) {
                    //l1rec->Lt[ipb] *= Fobar[ib] * l1rec->fsol * cos(l1rec->solz[ip]/OEL_RADEG) / OEL_PI;
                    l1rec->Lt[ipb] *= red_solar_irradiance[band] * cos(l1rec->solz[ip]/OEL_RADEG) / earth_sun_distance_correction / OEL_PI / 10.0;
                } else {
                    l1rec->Lt[ipb] /= 10.; // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the SWIR bands choose low gain if either high gain band saturates
        // band 2 low gain
        // band 3 high gain
        // band 5 low gain
        // band 6 high gain

        int saturated = 0;
        if((qual_SWIR[3][ip] & SATURATION_BIT) || (qual_SWIR[6][ip] & SATURATION_BIT)) {
            saturated = 1;
            l1rec->hilt[ip] = 1;
        }

        for(band=0; band<expected_num_SWIR_bands; band++) {

            if((band==2 || band==5) && !saturated)
                continue;
            else if((band==3 || band==6) && saturated)
                continue;

            if(Lt_SWIR[band][ip] == Lt_SWIRFillValue) {
                l1rec->Lt[ipb] = BAD_FLT;
            } else {
                l1rec->Lt[ipb] = Lt_SWIR[band][ip];
                if (use_rhot) {
                    //l1rec->Lt[ipb] *= Fobar[ib] * l1rec->fsol * cos(l1rec->solz[ip]/OEL_RADEG) / OEL_PI;
                    l1rec->Lt[ipb] *= SWIR_solar_irradiance[band] * cos(l1rec->solz[ip]/OEL_RADEG) / earth_sun_distance_correction / OEL_PI / 10.0;
                } else {
                    l1rec->Lt[ipb] /= 10.; // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }
    }

    // read ccd_scan_angles for current scan

    // read starting at current scan, starting at index 0
    start[0] = line;
    start[1] = 0;
    start[2] = 0;

    // read the line and num_pixels of data
    count[0] = 1;
    count[1] = num_pixels;
    count[2] = 1;

    // grab the id for ccd_scan_angle so that it can read from the file 
    status = nc_get_vara_float(navigationGrp, ccdScanAnglesId, start, count, polcorPrivateData->ccdScanAngle);
    check_err(status, __LINE__, __FILE__);

    return (LIFE_IS_GOOD);
}


/**
 * Close L1B file, GEO file, and free memory
 * @param file
 * @return 
 */
int closel1_oci(filehandle *file) {
    int status;

    status = nc_close(file->sd_id);
    check_err(status, __LINE__, __FILE__);

    // Free memory
    // From openl1_oci
    if (tmpShort) free(tmpShort);
    if (scan_time) free(scan_time);
    if (hamside) free(hamside);
    if (tilt) free(tilt);  
    if (scanQual) free(scanQual);
    if (Lt_blue) free2d_float(Lt_blue);
    if (Lt_red) free2d_float(Lt_red);
    if (Lt_SWIR) free2d_float(Lt_SWIR);
    if (qual_SWIR) free2d_uchar(qual_SWIR);
    if (blue_solar_irradiance) free(blue_solar_irradiance);
    if (red_solar_irradiance) free(red_solar_irradiance);
    if (SWIR_solar_irradiance) free(SWIR_solar_irradiance);

    return (LIFE_IS_GOOD);
}
