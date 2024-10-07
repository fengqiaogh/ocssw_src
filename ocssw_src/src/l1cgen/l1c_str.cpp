
//  l1c_str.cpp
//
//
//  Created by Martin Montes on 8/15/2022
#include "l1c_str.h"
#include <iostream>
#include <string>
#include "l1c_filehandle.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <stdint.h>
#include <netcdf.h>
#include <nc4utils.h>
#include "libnav.h"
#include <timeutils.h>
#include <genutils.h>
#include <allocate4d.h>
#include <allocate3d.h>
#include <allocate2d.h>
#include <l1.h>
#include <netcdf>

#include "l1c_latlongrid.h"

// static float *Fobar; // reflectance to radiance conversion factors
// static int extract_pixel_start = 0;

static short *tmpShort;
static int normalizedLt = 0;
static int use_rhot = 1;
// whole file stuff
/*
// ONLY FOR OCIS!!
static size_t expected_num_blue_bands = 120;  // 60;before
static size_t expected_num_red_bands = 120;   // 60;
static size_t expected_num_SWIR_bands = 9;
static size_t tot_num_bands = 239;  // 249-10
*/

//OCI
static size_t expected_num_blue_bands = 121;
static size_t expected_num_red_bands = 165;
static size_t expected_num_SWIR_bands = 9;
static size_t tot_num_bands = 295;  


// static size_t nviews;
// static size_t nbands;
static size_t num_scans, num_pixels;
static size_t num_blue_bands, num_red_bands, num_SWIR_bands, nband_view, npol_band_view;

static int ncid_L1B;

// scan line attributes
static double *scan_time;  // seconds of day
static int32_t scan_time_year, scan_time_month, scan_time_day;
// static int32_t scan_time_hour,scan_time_min;
// static double scan_time_secs;

// static uint_8 *scan_quality;

// geolocation data
static int geolocationGrp, lambdaGrp, blueGrp;
static int lonId, latId, senzId, senaId, solzId, solaId, bwId, rwId, swId, iwId, pwId, IntId, IpolId, vpId;
static float latFillValue =BAD_FLT;
static float lonFillValue = BAD_FLT;
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

static float I_FillValue = BAD_FLT;
static float Ipol_FillValue = BAD_FLT;

// Observation data
static int observationGrp;
static int Lt_blueId, Lt_redId, Lt_SWIRId;
// static float **Lt_blue;                //[num_blue_bands][num_pixels], This scan
// static float **Lt_red;                 //[num_red_bands][num_pixels], This scan
// static float **Lt_SWIR;                //[num_SWIR_bands][num_pixels], This scan

// static float ****Lt_all;

static float Lt_blueFillValue = BAD_FLT;
static float Lt_redFillValue = BAD_FLT;
static float Lt_SWIRFillValue = BAD_FLT;

// Navigation data
// static int navGrp;
// static int ovId,opId,attID;
// static float *orbv;//for each scan we have xpixels dim
// static float *orbp;
// static float *attang;
// static double *timeArr;

static int status, dimid;

using namespace std;

namespace l1c {

l1c_str::l1c_str()
    : att_ang{-999.0, -999.0, -999.0}, orb_pos{-999.0, -999.0, -999.0}, orb_vel{-999.0, -999.0, -999.0} {
    // global attributes
    npix = -1;
    iscan = 0;
    nscan = -1;
    nbands = -1;
    // scan line attributes--
    ev_mid_time = nullptr;
    scan_quality_flag = nullptr;
    spix = -1;
    epix = -1;
    dpix = -1;

    // structure--pointers to data arrays
    Lt = nullptr;       // dim depends between sensors--OCI is bands x pixels --
    Lt_blue = nullptr;  //[num_views][num_pol][num_blue_bands][num_pixels]
    Lt_red = nullptr;
    Lt_SWIR = nullptr;
    Lt_tot = nullptr;
    Fobar = nullptr;

    blue_lambdas = nullptr;
    red_lambdas = nullptr;
    SWIR_lambdas = nullptr;

    I = nullptr;
    I_polsample = nullptr;
    I_lambdas = nullptr;
    pol_lambdas = nullptr;
    viewport = nullptr;

    // sensor/sun geometry
    senz = nullptr;
    sena = nullptr;
    solz = nullptr;
    sola = nullptr;
    delphi = nullptr;
    scattang = nullptr;

    // geolocation--
    timepix = nullptr;
    senazpix = nullptr;
    latpix = nullptr;
    lonpix = nullptr;
    latpix2 = nullptr;
    lonpix2 = nullptr;
    latpix3 = nullptr;
    lonpix3 = nullptr;

    senazpix_3d = nullptr;
    latpix_3d = nullptr;
    lonpix_3d = nullptr;
    latpix2_3d = nullptr;  // lat2/lon are the sline +1
    lonpix2_3d = nullptr;

    latnad = nullptr;
    lonnad = nullptr;
    lonershift = nullptr;
    terr_height = nullptr;
    cloud_height = nullptr;

    l1cfile = nullptr;
    l1file = nullptr;
    // ancill info?
}

l1c_str::~l1c_str() {
}

int32_t l1c_str::openl1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file) {
    std::string str;
    const char *ptstr;
    // Open the netcdf4 input file
    str = l1cfile->l1b_name;
    ptstr = str.c_str();
    printf("Opening MISR L1B file\n");
    cout << ptstr << endl;

    if (openl1(l1file) != 0) {
        printf("-E- %s: Error opening for reading.\n", l1file->name);
        exit(FATAL_ERROR);
    }

    return 0;
}

int32_t l1c_str::openl1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    std::string str;
    const char *ptstr;
    // Open the netcdf4 input file
    str = l1cfile->l1b_name;
    ptstr = str.c_str();
    printf("Opening SPEXone L1B file\n");
    cout << str << endl;

    status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n", __FILE__, __LINE__, ptstr);
        return (1);
    }

    if (l1cfile->format == FT_SPEXONE) {
        // number of views
        status = nc_inq_dimid(ncid_L1B, "number_of_views", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_views.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &nviews);
        l1cstr->nviews = nviews;

        // number of scans
        status = nc_inq_dimid(ncid_L1B, "bins_along_track", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_scans.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);
        l1cstr->nscan = num_scans;

        // num_pixels
        status = nc_inq_dimid(ncid_L1B, "spatial_samples_per_image", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_pixels.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);
        l1cstr->npix = num_pixels;

        // number of bands
        status = nc_inq_dimid(ncid_L1B, "intensity_bands_per_view", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading intensity_bands_per_view.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &nband_view);

        status = nc_inq_dimid(ncid_L1B, "polarization_bands_per_view", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading polarization_bands_per_view.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &npol_band_view);
    }

    printf("L1B Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);

    // allocate all of the data
    l1cstr->senazpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->latpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->latpix2 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix2 = (float *)calloc(num_pixels, sizeof(float));

    l1cstr->I = allocate2d_float(num_pixels, nband_view);
    l1cstr->I_polsample = allocate2d_float(num_pixels, npol_band_view);
    l1cstr->I_lambdas = allocate2d_float(nviews, nband_view);
    l1cstr->pol_lambdas = allocate2d_float(nviews, npol_band_view);
    l1cstr->viewport = (uint8_t *)calloc(num_pixels, sizeof(uint8_t));

    // get start time
    size_t att_len;
    status = nc_inq_attlen(ncid_L1B, NC_GLOBAL, "time_coverage_start", &att_len);
    check_err(status, __LINE__, __FILE__);
    // allocate required space before retrieving values
    char *time_str = (char *)malloc(att_len + 1);  // + 1 for trailing null
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

    // read the orbit#
    int orbit_number;
    status = nc_get_att_int(ncid_L1B, NC_GLOBAL, "orbit_number", &orbit_number);
    check_err(status, __LINE__, __FILE__);

    // lambdas for each spectral range
    status = nc_inq_grp_ncid(ncid_L1B, "SENSOR_VIEW_BANDS", &lambdaGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "intensity_wavelengths", &iwId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "polarization_wavelengths", &pwId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "viewport_index", &vpId);
    check_err(status, __LINE__, __FILE__);

    status = nc_get_var_float(lambdaGrp, iwId, l1cstr->I_lambdas[0]);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, pwId, l1cstr->pol_lambdas[0]);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_ubyte(lambdaGrp, vpId, l1cstr->viewport);
    check_err(status, __LINE__, __FILE__);

    // Setup geofile pointers
    status = nc_inq_grp_ncid(ncid_L1B, "GEOLOCATION_DATA", &geolocationGrp);
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
    //  status = nc_get_att_float(geolocationGrp, senzId, "scale_factor", &senzScale);
    //    check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, senzId, "add_offset", &senzOffset);
    //    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "sensor_azimuth", &senaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, senaId, NULL, &senaFillValue);
    check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, senaId, "scale_factor", &senaScale);
    //    check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, senaId, "add_offset", &senaOffset);
    //   check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "solar_zenith", &solzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solzId, NULL, &solzFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solzId, "scale_factor", &solzScale);
    //   check_err(status, __LINE__, __FILE__);
    //   status = nc_get_att_float(geolocationGrp, solzId, "add_offset", &solzOffset);
    //    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "solar_azimuth", &solaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solaId, NULL, &solaFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, solaId, "scale_factor", &solaScale);
    //   check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, solaId, "add_offset", &solaOffset);
    //    check_err(status, __LINE__, __FILE__);

    // get IDs for the observations
    status = nc_inq_grp_ncid(ncid_L1B, "OBSERVATION_DATA", &observationGrp);
    check_err(status, __LINE__, __FILE__);

    // Get varids for each of the Lt_*
    status = nc_inq_varid(observationGrp, "I", &IntId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(observationGrp, IntId, NULL, &I_FillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(observationGrp, "I_polsample", &IpolId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(observationGrp, IpolId, NULL, &Ipol_FillValue);
    check_err(status, __LINE__, __FILE__);

    // l1cfile assig
    l1cfile->sd_id = ncid_L1B;
    l1cfile->nband_view = nband_view;
    l1cfile->npol_band_view = npol_band_view;
    l1cfile->n_views = nviews;

    cout << "nband_view =" << nband_view << endl;
    cout << "npol_band_view = " << npol_band_view << endl;

    l1cfile->npix = num_pixels;
    l1cfile->nadpix = (num_pixels - 1) / 2;  // nadir pixel index
    l1cfile->nscan = num_scans;
    l1cfile->terrain_corrected = 1;  // presumed.
    l1cfile->orbit_number = orbit_number;

    cout << "n_views..." << l1cfile->n_views << "intensity bands.per view..." << nband_view
         << "polarization bands.per view." << npol_band_view << endl;

    return 0;
}

int32_t l1c_str::openl1b_oci_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    std::string str;
    const char *ptstr;

    // Open the netcdf4 input file
    str = l1cfile->l1b_name;
    ptstr = str.c_str();
    printf("Opening OCI L1B file\n");
    cout << str << endl;

    status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n", __FILE__, __LINE__, ptstr);
        return (1);
    }

    // num_scans
    if (l1cfile->format == FT_OCIL1B) {
        status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_scans.\n");
            exit(EXIT_FAILURE);
        }
        cout << "reading OCI" << endl;
        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);
        l1cstr->nscan = num_scans;

        // num_pixels
        status = nc_inq_dimid(ncid_L1B, "ccd_pixels", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_pixels.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);
        l1cstr->npix = num_pixels;

        // num_blue_bands
        status = nc_inq_dimid(ncid_L1B, "blue_bands", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_blue_bands.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_blue_bands);
        if (num_blue_bands < expected_num_blue_bands) {
            fprintf(stderr, "-E- blue bands in file are less than expected!, expecting %d, found %d.\n",
                    (int)expected_num_blue_bands, (int)num_blue_bands);
//            exit(EXIT_FAILURE);
        }

        // num_red_bands
        status = nc_inq_dimid(ncid_L1B, "red_bands", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_red_bands.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_red_bands);
        if (num_red_bands < expected_num_red_bands) {
            fprintf(stderr, "-E- red bands in file are less than expected!, expecting %d, found %d.\n",
                    (int)expected_num_red_bands, (int)num_red_bands);
//            exit(EXIT_FAILURE);
        }

        // num_SWIR_bands
        status = nc_inq_dimid(ncid_L1B, "SWIR_bands", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_SWIR_bands.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_SWIR_bands);
        if (num_SWIR_bands < expected_num_SWIR_bands) {
            fprintf(stderr, "-E- SWIR bands in file are less than expected!, expecting %d, found %d.\n",
                    (int)expected_num_SWIR_bands, (int)num_SWIR_bands);
//            exit(EXIT_FAILURE);
        }
    }

    printf("L1B Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);

    // allocate all of the data
    // line by line!!!
    tmpShort = (short *)calloc(num_pixels, sizeof(short));
    l1cstr->timepix = (double *)calloc(num_scans, sizeof(double));
    l1cstr->tilt = (float *)calloc(num_scans, sizeof(float));
    l1cstr->latpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->latpix2 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix2 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->latpix3 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix3 = (float *)calloc(num_pixels, sizeof(float));

    l1cstr->Lt_blue = allocate2d_float(num_blue_bands, num_pixels);
    l1cstr->Lt_red = allocate2d_float(num_red_bands, num_pixels);
    l1cstr->Lt_SWIR = allocate2d_float(num_SWIR_bands, num_pixels);
    l1cstr->Lt_tot = (float *)calloc(num_pixels * tot_num_bands, sizeof(float));

    l1cstr->Fobar = (float *)calloc(tot_num_bands, sizeof(float));
    l1cstr->solz = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->senz = (float *)calloc(num_pixels, sizeof(float));

    l1cstr->blue_lambdas = (float *)calloc(num_blue_bands, sizeof(float));
    l1cstr->red_lambdas = (float *)calloc(num_red_bands, sizeof(float));
    l1cstr->SWIR_lambdas = (float *)calloc(num_SWIR_bands, sizeof(float));

    scan_time = (double *)malloc(num_scans * sizeof(double));

    // Get group id from L1B file for GROUP scan_line_attributes.
    int groupid;
    int varid;

    if ((nc_inq_grp_ncid(ncid_L1B, "navigation_data", &groupid)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding navigation_data.\n");
        exit(EXIT_FAILURE);
    }
    status = nc_inq_varid(groupid, "tilt", &varid);
    float tiltFillValue = BAD_FLT;
    if (status == NC_NOERR) {
        status = nc_inq_var_fill(groupid, varid, NULL, &tiltFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(groupid, varid, l1cstr->tilt);
        check_err(status, __LINE__, __FILE__);
    }

    if ((nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &groupid)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding scan_line_attributes.\n");
        exit(EXIT_FAILURE);
    }

    double scan_timeFillValue = BAD_FLT;
    status = nc_inq_varid(groupid, "time", &varid);
    if (status == NC_NOERR) {
        status = nc_inq_var_fill(groupid, varid, NULL, &scan_timeFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(groupid, varid, l1cstr->timepix);
        check_err(status, __LINE__, __FILE__);
        /*       status = nc_get_att_int(groupid, varid, "year", &scan_time_year);
               check_err(status, __LINE__, __FILE__);
               status = nc_get_att_int(groupid, varid, "month", &scan_time_month);
               check_err(status, __LINE__, __FILE__);
               status = nc_get_att_int(groupid, varid, "day", &scan_time_day);
               check_err(status, __LINE__, __FILE__);
        */
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
        char *time_str = (char *)malloc(att_len + 1);  // + 1 for trailing null

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

    for (size_t i = 0; i < num_scans; i++) {
        if (scan_time[i] == scan_timeFillValue)
            scan_time[i] = BAD_FLT;
    }

    // read the orbit#
    /*    int orbit_number;
        status = nc_get_att_int(ncid_L1B, NC_GLOBAL, "orbit_number", &orbit_number);
        check_err(status, __LINE__, __FILE__);
    */
    // lambdas for each spectral range
    status = nc_inq_grp_ncid(ncid_L1B, "sensor_band_parameters", &lambdaGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "blue_wavelength", &bwId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "red_wavelength", &rwId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "SWIR_wavelength", &swId);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, bwId, l1cstr->blue_lambdas);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, rwId, l1cstr->red_lambdas);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, swId, l1cstr->SWIR_lambdas);
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
    /*    status = nc_get_att_float(geolocationGrp, senzId, "scale_factor", &senzScale);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_att_float(geolocationGrp, senzId, "add_offset", &senzOffset);
        check_err(status, __LINE__, __FILE__);
    */
    status = nc_inq_varid(geolocationGrp, "sensor_azimuth", &senaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, senaId, NULL, &senaFillValue);
    check_err(status, __LINE__, __FILE__);
    /*    status = nc_get_att_float(geolocationGrp, senaId, "scale_factor", &senaScale);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_att_float(geolocationGrp, senaId, "add_offset", &senaOffset);
        check_err(status, __LINE__, __FILE__);
    */
    status = nc_inq_varid(geolocationGrp, "solar_zenith", &solzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solzId, NULL, &solzFillValue);
    check_err(status, __LINE__, __FILE__);
    /*    status = nc_get_att_float(geolocationGrp, solzId, "scale_factor", &solzScale);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_att_float(geolocationGrp, solzId, "add_offset", &solzOffset);
        check_err(status, __LINE__, __FILE__);
    */
    status = nc_inq_varid(geolocationGrp, "solar_azimuth", &solaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, solaId, NULL, &solaFillValue);
    check_err(status, __LINE__, __FILE__);
    /*    status = nc_get_att_float(geolocationGrp, solaId, "scale_factor", &solaScale);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_att_float(geolocationGrp, solaId, "add_offset", &solaOffset);
        check_err(status, __LINE__, __FILE__);
    */

    // get IDs for the observations
    status = nc_inq_grp_ncid(ncid_L1B, "observation_data", &observationGrp);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(observationGrp, "Lt_blue", &Lt_blueId);

    use_rhot = 0;
    // Get varids for each of the Lt_*
    if (status == NC_NOERR) {
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
    } else {
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

    l1cstr->nbands = tot_num_bands;
    // l1cfile assig
    l1cfile->sd_id = ncid_L1B;
    l1cfile->nband_blue = num_blue_bands;
    l1cfile->nband_red = num_red_bands;
    l1cfile->nband_swir = num_SWIR_bands;
    l1cfile->nbands = num_blue_bands + num_red_bands + num_SWIR_bands;
    nviews = 2;
    l1cfile->n_views = nviews;

    cout << "file->nbands_b =" << l1cfile->nband_blue << endl;
    cout << "file->nbands_r = " << l1cfile->nband_red << endl;
    cout << "file->nbands_swir = " << l1cfile->nband_swir << endl;
    /*
        if(l1cfile->format==FT_OCIL1B){
            expected_num_blue_bands=120;//actually 121,165,9
            expected_num_red_bands=168;
            expected_num_SWIR_bands=9;
            tot_num_bands=expected_num_blue_bands+expected_num_red_bands+expected_num_SWIR_bands;
        }
    */
    //    cout<<"expected blue.."<<expected_num_blue_bands<<"expected
    //    red.."<<expected_num_red_bands<<"expected SWIR.."<<expected_num_SWIR_bands<<"tot
    //    bands.."<<tot_num_bands<<endl;

    l1cfile->npix = num_pixels;
    l1cfile->nadpix = (num_pixels - 1) / 2;  // nadir pixel index
    l1cfile->nscan = num_scans;
    l1cfile->ndets = 1;
    l1cfile->terrain_corrected = 1;  // presumed.
    //    l1cfile->orbit_number = orbit_number;

    return 0;
}

int32_t l1c_str::openl1b_ocis_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    std::string str;
    const char *ptstr;

    // Open the netcdf4 input file
    str = l1cfile->l1b_name;
    ptstr = str.c_str();
    printf("Opening OCIS L1B file\n");
    cout << str << endl;

    status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n", __FILE__, __LINE__, ptstr);
        return (1);
    }

    // num_scans
    if (l1cfile->format == FT_OCIS) {
        status = nc_inq_dimid(ncid_L1B, "number_of_scans", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_scans.\n");
            exit(EXIT_FAILURE);
        }
        cout << "reading OCIS" << endl;
        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);
        l1cstr->nscan = num_scans;

        // num_pixels
        status = nc_inq_dimid(ncid_L1B, "ccd_pixels", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_pixels.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);
        l1cstr->npix = num_pixels;

        // num_blue_bands
        status = nc_inq_dimid(ncid_L1B, "blue_bands", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_blue_bands.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_blue_bands);
        if (num_blue_bands < expected_num_blue_bands) {
            fprintf(stderr, "-E- blue bands in file is less than expected!, expecting %d, found %d.\n",
                    (int)expected_num_blue_bands, (int)num_blue_bands);
       //     exit(EXIT_FAILURE);
        }

        // num_red_bands
        status = nc_inq_dimid(ncid_L1B, "red_bands", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_red_bands.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_red_bands);
        if (num_red_bands < expected_num_red_bands) {
            fprintf(stderr, "-E- red bands in file is less than expected!, expecting %d, found %d.\n",
                    (int)expected_num_red_bands, (int)num_red_bands);
  //          exit(EXIT_FAILURE);
        }

        // num_SWIR_bands
        status = nc_inq_dimid(ncid_L1B, "SWIR_bands", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_SWIR_bands.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_SWIR_bands);
        if (num_SWIR_bands < expected_num_SWIR_bands) {
            fprintf(stderr, "-E- SWIR bands in file is less than expected!, expecting %d, found %d.\n",
                    (int)expected_num_SWIR_bands, (int)num_SWIR_bands);
//            exit(EXIT_FAILURE);
        }
    }

    printf("L1B Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);

    // allocate all of the data
    // line by line!!!
    tmpShort = (short *)calloc(num_pixels, sizeof(short));
    l1cstr->timepix = (double *)calloc(num_scans, sizeof(double));
    l1cstr->tilt = (float *)calloc(num_scans, sizeof(float));
    l1cstr->latpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->latpix2 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix2 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->latpix3 = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->lonpix3 = (float *)calloc(num_pixels, sizeof(float));

    l1cstr->Lt_blue = allocate2d_float(num_blue_bands, num_pixels);
    l1cstr->Lt_red = allocate2d_float(num_red_bands, num_pixels);
    l1cstr->Lt_SWIR = allocate2d_float(num_SWIR_bands, num_pixels);
    l1cstr->Lt_tot = (float *)calloc(num_pixels * tot_num_bands, sizeof(float));

    l1cstr->Fobar = (float *)calloc(tot_num_bands, sizeof(float));
    l1cstr->solz = (float *)calloc(num_pixels, sizeof(float));
    l1cstr->senz = (float *)calloc(num_pixels, sizeof(float));

    l1cstr->blue_lambdas = (float *)calloc(num_blue_bands, sizeof(float));
    l1cstr->red_lambdas = (float *)calloc(num_red_bands, sizeof(float));
    l1cstr->SWIR_lambdas = (float *)calloc(num_SWIR_bands, sizeof(float));

    scan_time = (double *)malloc(num_scans * sizeof(double));

    // Get group id from L1B file for GROUP scan_line_attributes.
    int groupid;
    int varid;

    if ((nc_inq_grp_ncid(ncid_L1B, "navigation_data", &groupid)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding navigation_data.\n");
        exit(EXIT_FAILURE);
    }
    status = nc_inq_varid(groupid, "tilt_angle", &varid);
    float tiltFillValue = BAD_FLT;
    if (status == NC_NOERR) {
        status = nc_inq_var_fill(groupid, varid, NULL, &tiltFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_float(groupid, varid, l1cstr->tilt);
        check_err(status, __LINE__, __FILE__);
    }

    if ((nc_inq_grp_ncid(ncid_L1B, "scan_line_attributes", &groupid)) == NC_NOERR) {
    } else {
        fprintf(stderr, "-E- Error finding scan_line_attributes.\n");
        exit(EXIT_FAILURE);
    }

    double scan_timeFillValue = BAD_FLT;
    status = nc_inq_varid(groupid, "time", &varid);
    if (status == NC_NOERR) {
        status = nc_inq_var_fill(groupid, varid, NULL, &scan_timeFillValue);
        check_err(status, __LINE__, __FILE__);
        status = nc_get_var_double(groupid, varid, l1cstr->timepix);
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
        char *time_str = (char *)malloc(att_len + 1);  // + 1 for trailing null

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

    for (size_t i = 0; i < num_scans; i++) {
        if (scan_time[i] == scan_timeFillValue)
            scan_time[i] = BAD_FLT;
    }

    // read the orbit#
    int orbit_number;
    status = nc_get_att_int(ncid_L1B, NC_GLOBAL, "orbit_number", &orbit_number);
    check_err(status, __LINE__, __FILE__);

    // lambdas for each spectral range
    status = nc_inq_grp_ncid(ncid_L1B, "sensor_band_parameters", &lambdaGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "blue_wavelength", &bwId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "red_wavelength", &rwId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(lambdaGrp, "SWIR_wavelength", &swId);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, bwId, l1cstr->blue_lambdas);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, rwId, l1cstr->red_lambdas);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_var_float(lambdaGrp, swId, l1cstr->SWIR_lambdas);
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

    status = nc_inq_varid(observationGrp, "Lt_blue", &Lt_blueId);

    use_rhot = 0;
    // Get varids for each of the Lt_*
    if (status == NC_NOERR) {
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
    } else {
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

    l1cstr->nbands = tot_num_bands;
    // l1cfile assig
    l1cfile->sd_id = ncid_L1B;
    l1cfile->nband_blue = num_blue_bands;
    l1cfile->nband_red = num_red_bands;
    l1cfile->nband_swir = num_SWIR_bands;
    l1cfile->nbands = num_blue_bands + num_red_bands + num_SWIR_bands;
    nviews = 2;
    l1cfile->n_views = nviews;

    cout << "file->nbands_b =" << l1cfile->nband_blue << endl;
    cout << "file->nbands_r = " << l1cfile->nband_red << endl;
    cout << "file->nbands_swir = " << l1cfile->nband_swir << endl;
    /*
        if(l1cfile->format==FT_OCIL1B){
            expected_num_blue_bands=120;//actually 121,165,9
            expected_num_red_bands=168;
            expected_num_SWIR_bands=9;
            tot_num_bands=expected_num_blue_bands+expected_num_red_bands+expected_num_SWIR_bands;
        }
    */
    //   cout<<"expected blue.."<<expected_num_blue_bands<<"expected red.."<<expected_num_red_bands<<"expected
    //   SWIR.."<<expected_num_SWIR_bands<<"tot bands.."<<tot_num_bands<<endl;

    l1cfile->npix = num_pixels;
    l1cfile->nadpix = (num_pixels - 1) / 2;  // nadir pixel index
    l1cfile->nscan = num_scans;
    l1cfile->ndets = 1;
    l1cfile->terrain_corrected = 1;  // presumed.
    //    l1cfile->orbit_number = orbit_number;

    return 0;
}

int32_t l1c_str::openl1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    std::string str;
    const char *ptstr;

    // Open the netcdf4 input file
    str = l1cfile->l1b_name;
    ptstr = str.c_str();
    printf("Opening HARP2 L1B file\n");
    cout << str << endl;

    status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n", __FILE__, __LINE__, ptstr);
        return (1);
    }

    if (l1cfile->format == FT_HARP2) {
        // open blue group--
        status = nc_inq_grp_ncid(ncid_L1B, "blue", &blueGrp);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading group blue.\n");
            exit(EXIT_FAILURE);
        }

        // number of views
        status = nc_inq_dimid(blueGrp, "Views", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_views.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &nviews);
        l1cstr->nviews = nviews;

        // number of scans
        status = nc_inq_dimid(ncid_L1B, "Swath_Lines", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading number_of_lines.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_scans);
        l1cstr->nscan = num_scans;

        // num_pixels
        status = nc_inq_dimid(ncid_L1B, "Swath_Pixels", &dimid);
        if (status != NC_NOERR) {
            fprintf(stderr, "-E- Error reading num_pixels.\n");
            exit(EXIT_FAILURE);
        }
        nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);
        l1cstr->npix = num_pixels;

        /*       //number of bands
                status = nc_inq_dimid(ncid_L1B, "intensity_bands_per_view", &dimid);
                if (status != NC_NOERR) {
                   fprintf(stderr, "-E- Error reading intensity_bands_per_view.\n");
                   exit(EXIT_FAILURE);
                  }
               nc_inq_dimlen(ncid_L1B, dimid, &nband_view);


                status = nc_inq_dimid(ncid_L1B, "polarization_bands_per_view", &dimid);
                if (status != NC_NOERR) {
                   fprintf(stderr, "-E- Error reading polarization_bands_per_view.\n");
                   exit(EXIT_FAILURE);
                  }
               nc_inq_dimlen(ncid_L1B, dimid, &npol_band_view);
        */

        // blue group
        nband_view = 4;
        npol_band_view = 4;
    }

    printf("L1B Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);

    // allocate all of the data
    l1cstr->senazpix_3d = allocate2d_float(nviews, num_pixels);
    l1cstr->latpix_3d = allocate2d_float(nviews, num_pixels);
    l1cstr->lonpix_3d = allocate2d_float(nviews, num_pixels);
    l1cstr->latpix2_3d = allocate2d_float(nviews, num_pixels);
    l1cstr->lonpix2_3d = allocate2d_float(nviews, num_pixels);

    l1cstr->I = allocate2d_float(nviews, num_pixels);
    l1cstr->I_polsample = allocate2d_float(nviews, num_pixels);

    l1cstr->I_lambdas = allocate2d_float(nviews, nband_view);
    l1cstr->pol_lambdas = allocate2d_float(nviews, npol_band_view);

    // only blue group
    for (size_t i = 0; i < nviews; i++) {
        l1cstr->I_lambdas[i][0] = 441.9;
        l1cstr->I_lambdas[i][1] = 549.8;
        l1cstr->I_lambdas[i][2] = 669.4;
        l1cstr->I_lambdas[i][3] = 867.8;

        l1cstr->pol_lambdas[i][0] = 441.9;
        l1cstr->pol_lambdas[i][1] = 549.8;
        l1cstr->pol_lambdas[i][2] = 669.4;
        l1cstr->pol_lambdas[i][3] = 867.8;
    }

    // Get group id from L1B file for GROUP scan_line_attributes.

    // get start time
    size_t att_len;
    status = nc_inq_attlen(ncid_L1B, NC_GLOBAL, "time_coverage_start", &att_len);
    check_err(status, __LINE__, __FILE__);

    // allocate required space before retrieving values
    char *time_str = (char *)malloc(att_len + 1);  // + 1 for trailing null

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

    /*
       for(size_t i=0; i<num_scans; i++) {
          if(scan_time[i] == scan_timeFillValue)
              scan_time[i] = BAD_FLT;
      }
  */

    // read the orbit#
    int orbit_number;
    status = nc_get_att_int(ncid_L1B, NC_GLOBAL, "orbit_number", &orbit_number);
    check_err(status, __LINE__, __FILE__);

    // lambdas for each spectral range
    /*
       status = nc_inq_grp_ncid(ncid_L1B, "SENSOR_VIEW_BANDS", &lambdaGrp);
       check_err(status, __LINE__, __FILE__);
       status = nc_inq_varid(lambdaGrp, "intensity_wavelengths", &iwId);
       check_err(status, __LINE__, __FILE__);
       status = nc_inq_varid(lambdaGrp, "polarization_wavelengths", &pwId);
       check_err(status, __LINE__, __FILE__);




       status = nc_get_var_float(lambdaGrp, iwId,l1cstr->I_lambdas[0]);
           check_err(status, __LINE__, __FILE__);
       status = nc_get_var_float(lambdaGrp, pwId,l1cstr->pol_lambdas[0]);
           check_err(status, __LINE__, __FILE__);
   */

    // Setup geofile pointers
    status = nc_inq_varid(blueGrp, "Longitude", &lonId);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_var_fill(blueGrp, lonId, NULL, &lonFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(blueGrp, "Latitude", &latId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(blueGrp, latId, NULL, &latFillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(blueGrp, "View_Zenith", &senzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(blueGrp, senzId, NULL, &senzFillValue);
    check_err(status, __LINE__, __FILE__);
    //  status = nc_get_att_float(geolocationGrp, senzId, "scale_factor", &senzScale);
    //    check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, senzId, "add_offset", &senzOffset);
    //    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(blueGrp, "View_Azimuth", &senaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(blueGrp, senaId, NULL, &senaFillValue);
    check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, senaId, "scale_factor", &senaScale);
    //    check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, senaId, "add_offset", &senaOffset);
    //   check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(blueGrp, "Solar_Zenith", &solzId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(blueGrp, solzId, NULL, &solzFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(blueGrp, solzId, "scale_factor", &solzScale);
    //   check_err(status, __LINE__, __FILE__);
    //   status = nc_get_att_float(geolocationGrp, solzId, "add_offset", &solzOffset);
    //    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(blueGrp, "Solar_Azimuth", &solaId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(blueGrp, solaId, NULL, &solaFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(blueGrp, solaId, "scale_factor", &solaScale);
    //   check_err(status, __LINE__, __FILE__);
    //    status = nc_get_att_float(geolocationGrp, solaId, "add_offset", &solaOffset);
    //    check_err(status, __LINE__, __FILE__);

    // get IDs for the observations
    // Get varids for each of the Lt_*
    status = nc_inq_varid(blueGrp, "I", &IntId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(blueGrp, IntId, NULL, &I_FillValue);
    check_err(status, __LINE__, __FILE__);
    /*
        status = nc_inq_varid(blueGrp, "I_polsample", &IpolId);
        check_err(status, __LINE__, __FILE__);
        status = nc_inq_var_fill(blueGrp, IpolId, NULL, &Ipol_FillValue);
        check_err(status, __LINE__, __FILE__);
    */

    // l1cfile assig
    l1cfile->sd_id = ncid_L1B;
    l1cfile->nband_view = nband_view;
    l1cfile->npol_band_view = npol_band_view;
    l1cfile->n_views = nviews;

    cout << "nband_view =" << nband_view << endl;
    cout << "npol_band_view = " << npol_band_view << endl;

    l1cfile->npix = num_pixels;
    l1cfile->nadpix = (num_pixels - 1) / 2;  // nadir pixel index
    l1cfile->nscan = num_scans;
    l1cfile->terrain_corrected = 1;  // presumed.
    l1cfile->orbit_number = orbit_number;

    cout << "n_views..." << l1cfile->n_views << "intensity bands.per view..." << nband_view
         << "polarization bands.per view." << npol_band_view << endl;

    return 0;
}

//--------------------------------------------------------------------------------------
int32_t l1c_str::readl1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t sline) {
    size_t oneline;
    size_t start[] = {0, 0}, start2[] = {0, 0}, start3[] = {0, 0, 0};
    size_t count[] = {1, 1}, count2[] = {1, 1, 1};
    string str;

    str = l1cfile->l1b_name;

    printf("Reading SPEXone L1B file\n");
    cout << str << endl;

    // assigning mem for Lt of different bands **********************

    // GEOLOCATION
    l1cstr->iscan = sline;

    start[0] = sline;
    start[1] = 0;
    count[0] = 1;
    count[1] = num_pixels;  // 1 line at a time

    status = nc_get_vara_float(
        geolocationGrp, senaId, start, count,
        l1cstr->senazpix);  // sensor azimuth -1800 to 1800 degrees, so must be divdied by 100
    status = nc_get_vara_float(geolocationGrp, latId, start, count, l1cstr->latpix);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start, count, l1cstr->lonpix);
    check_err(status, __LINE__, __FILE__);

    oneline = sline;

    if (oneline < num_scans - 1)
        start2[0] = sline + 1;
    else
        start2[0] = sline;

    start2[1] = 0;

    status = nc_get_vara_float(geolocationGrp, latId, start2, count, l1cstr->latpix2);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start2, count, l1cstr->lonpix2);
    check_err(status, __LINE__, __FILE__);

    // RADIANCES
    start3[0] = sline;
    start3[1] = 0;
    start3[2] = 0;
    count2[0] = 1;
    count2[1] = num_pixels;  // 1 line at a time
    count2[2] = nband_view;

    status = nc_get_vara_float(observationGrp, IntId, start3, count2,
                               l1cstr->I[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    count2[2] = npol_band_view;
    status = nc_get_vara_float(observationGrp, IpolId, start3, count2, l1cstr->I_polsample[0]);
    check_err(status, __LINE__, __FILE__);

    return 0;
}

int32_t l1c_str::writel1c_ocis(l1c_str *l1cstr, bin_str *binl1c, netCDF::NcFile *nc_output, float **Ltfrac,
                               float **areafrac, short **obs_view, int band_ix, int view_ix) {
    string str;
    short minval, maxval;
    float minval2, maxval2;

    netCDF::NcDim yd = nc_output->getDim("bins_along_track");
    int ybins = yd.getSize();
    netCDF::NcDim xd = nc_output->getDim("bins_across_track");
    int xbins = xd.getSize();

    cout << "saving full L1C granule -- aw binning.."
         << "ybins.." << ybins << "xbins.." << xbins << "band #:  " << band_ix + 1
         << "view #: " << view_ix + 1 << endl;

    std::vector<size_t> start3;
    start3.push_back(0);
    start3.push_back(0);
    start3.push_back(view_ix);

    std::vector<size_t> count3;
    count3.push_back(ybins);
    count3.push_back(xbins);
    count3.push_back(1);

    short **obs_clean = allocate2d_short(ybins, xbins);

    netCDF::NcGroup od_grp = nc_output->getGroup("observation_data");
    netCDF::NcVar v1 = od_grp.getVar("obs_per_view");
    netCDF::NcVarAtt a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval);

    for (int i = 0; i < ybins; i++) {
        for (int j = 0; j < xbins; j++) {
            if (obs_view[i][j] > maxval || obs_view[i][j] < minval) {
                obs_view[i][j] = binl1c->fillval1;  // fillvalue
            }
            obs_clean[i][j] = obs_view[i][j];
        }
    }
    v1.putVar(start3, count3, &obs_clean[0][0]);
    free (obs_clean);

    v1 = od_grp.getVar("I");
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval2);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval2);

    std::vector<size_t> start4;
    start4.push_back(0);
    start4.push_back(0);
    start4.push_back(view_ix);
    start4.push_back(band_ix);

    std::vector<size_t> count4;
    count4.push_back(ybins);
    count4.push_back(xbins);
    count4.push_back(1);
    count4.push_back(1);

    float **Lt = allocate2d_float(ybins, xbins);

    for (int i = 0; i < ybins; i++) {
        for (int j = 0; j < xbins; j++) {
            if (areafrac[i][j] > 0) {
                Lt[i][j] = Ltfrac[i][j] / sqrt(areafrac[i][j]);

                // check min/max values
                if (Lt[i][j] < minval2 || Lt[i][j] > maxval2)
                    Lt[i][j] = binl1c->fillval2;
            } else {
                Lt[i][j] = binl1c->fillval2;  // fillvalue
            }
        }
    }

    v1.putVar(start4, count4, &Lt[0][0]);

    free (Lt);

    return 0;
}

int32_t l1c_str::readl1b_ocis_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t sline) {
    size_t oneline;
    size_t start[] = {0, 0}, start2[] = {0, 0}, start4[] = {0, 0, 0};
    size_t count[] = {1, 1}, count2[] = {1, 1, 1};
    string str;

    str = l1cfile->l1b_name;

    printf("Reading OCIS L1B file\n");
    cout << str << endl;

    // assigning mem for Lt of different bands **********************

    // GEOLOCATION
    l1cstr->iscan = sline;

    start[0] = sline;
    start[1] = 0;
    count[0] = 1;
    count[1] = num_pixels;  // 1 line at a time

    status = nc_get_vara_float(geolocationGrp, solzId, start, count, l1cstr->solz);
    status = nc_get_vara_float(geolocationGrp, senzId, start, count, l1cstr->senz);
    status = nc_get_vara_float(geolocationGrp, latId, start, count, l1cstr->latpix);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start, count, l1cstr->lonpix);
    check_err(status, __LINE__, __FILE__);

    oneline = sline;

    if (oneline < num_scans - 1)
        start2[0] = sline + 1;
    else
        start2[0] = sline;

    start2[1] = 0;

    status = nc_get_vara_float(geolocationGrp, latId, start2, count, l1cstr->latpix2);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start2, count, l1cstr->lonpix2);
    check_err(status, __LINE__, __FILE__);

    // RADIANCES---Don'way for indexing
    start4[0] = 0;
    start4[1] = sline;
    start4[2] = 0;
    count2[0] = num_blue_bands;
    count2[1] = 1;  // 1 line at a time
    count2[2] = num_pixels;

    status = nc_get_vara_float(observationGrp, Lt_blueId, start4, count2,
                               l1cstr->Lt_blue[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    count2[0] = num_red_bands;
    status = nc_get_vara_float(observationGrp, Lt_redId, start4, count2, l1cstr->Lt_red[0]);
    check_err(status, __LINE__, __FILE__);

    count2[0] = num_SWIR_bands;
    status = nc_get_vara_float(observationGrp, Lt_SWIRId, start4, count2, l1cstr->Lt_SWIR[0]);
    check_err(status, __LINE__, __FILE__);

    double scantime = ymds2unix((int16_t)scan_time_year, (int16_t)scan_time_month, (int16_t)scan_time_day,
                                scan_time[sline]);
    int16_t syear, sday;
    double secs;
    unix2yds(scantime, &syear, &sday, &secs);

    int32_t yr = syear;
    int32_t dy = sday;
    int32_t msec = (int32_t)(secs * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);

    float fsol = pow(1.0 / esdist, 2);

    for (size_t ip = 0; ip < num_pixels; ip++) {
        size_t band;
        int ib = 0;
        int ipb = ip * tot_num_bands;

        // load up the blue bands skip the last two
        if(num_blue_bands<expected_num_blue_bands){
            expected_num_blue_bands=num_blue_bands;
        }
        else{
            expected_num_blue_bands=expected_num_blue_bands-2;
        }

        for (band = 0; band < expected_num_blue_bands; band++) {
            if (Lt_blue[band][ip] == Lt_blueFillValue) {
                l1cstr->Lt_tot[ipb] = 0.001;  // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1cstr->Lt_tot[ipb] = Lt_blue[band][ip];
                if (normalizedLt) {
                    l1cstr->Lt_tot[ipb] *= 100;
                } else if (use_rhot) {
                    //    l1cstr->Lt_tot[ipb] *= l1cstr->Fobar[ib] * fsol * cos(l1cstr->solz[ip]/RADEG) /
                    //    M_PI;

                } else {
                    l1cstr->Lt_tot[ipb] /= 10.;  // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the red bands skipping the first two and the last 4
        if(num_red_bands<expected_num_red_bands){
            expected_num_red_bands=num_red_bands;
        }
        else{
            expected_num_red_bands=expected_num_red_bands-4;
        }

        for (band = 2; band < expected_num_red_bands; band++) {
            if (Lt_red[band][ip] == Lt_redFillValue) {
                l1cstr->Lt_tot[ipb] = 0.001;  // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1cstr->Lt_tot[ipb] = Lt_red[band][ip];
                if (normalizedLt) {
                    l1cstr->Lt_tot[ipb] *= 100;
                } else if (use_rhot) {
                    l1cstr->Lt_tot[ipb] *= l1cstr->Fobar[ib] * fsol * cos(l1cstr->solz[ip] / RADEG) / M_PI;
                } else {
                    l1cstr->Lt_tot[ipb] /= 10.;  // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the SWIR bands, skip band 3 and 6, hi/low gain wavelengths
        if(num_SWIR_bands<expected_num_SWIR_bands){
            expected_num_SWIR_bands=num_SWIR_bands;
        }

        for (band = 0; band < expected_num_SWIR_bands; band++) {
            if (band == 3 || band == 6)
                continue;
            if (Lt_SWIR[band][ip] == Lt_SWIRFillValue) {
                l1cstr->Lt_tot[ipb] = 0.001;  // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1cstr->Lt_tot[ipb] = Lt_SWIR[band][ip];
                if (normalizedLt) {
                    l1cstr->Lt_tot[ipb] *= 100;
                } else if (use_rhot) {
                    l1cstr->Lt_tot[ipb] *= l1cstr->Fobar[ib] * fsol * cos(l1cstr->solz[ip] / RADEG) / M_PI;
                } else {
                    l1cstr->Lt_tot[ipb] /= 10.;  // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }
    }

    return 0;
}

int32_t l1c_str::readl1b_ocis_3lines(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t sline) {
    size_t oneline;
    size_t start[] = {0, 0}, start2[] = {0, 0}, start3[] = {0, 0}, start4[] = {0, 0, 0};
    size_t count[] = {1, 1}, count2[] = {1, 1, 1};
    string str;

    str = l1cfile->l1b_name;

    printf("Reading OCIS L1B file\n");


    // assigning mem for Lt of different bands **********************

    // GEOLOCATION
    l1cstr->iscan = sline;

    start[0] = sline;
    start[1] = 0;
    count[0] = 1;
    count[1] = num_pixels;  // 1 line at a time

    status = nc_get_vara_float(geolocationGrp, solzId, start, count, l1cstr->solz);
    status = nc_get_vara_float(geolocationGrp, senzId, start, count, l1cstr->senz);
    status = nc_get_vara_float(geolocationGrp, latId, start, count, l1cstr->latpix);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start, count, l1cstr->lonpix);
    check_err(status, __LINE__, __FILE__);

    oneline = sline;

    if (oneline < num_scans - 1)
        start2[0] = sline + 1;
    else
        start2[0] = sline;

    start2[1] = 0;

    status = nc_get_vara_float(geolocationGrp, latId, start2, count, l1cstr->latpix2);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start2, count, l1cstr->lonpix2);
    check_err(status, __LINE__, __FILE__);

    if (oneline < num_scans - 2)
        start3[0] = sline + 2;
    else
        start3[0] = sline;

    start3[1] = 0;

    status = nc_get_vara_float(geolocationGrp, latId, start3, count, l1cstr->latpix3);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start3, count, l1cstr->lonpix3);
    check_err(status, __LINE__, __FILE__);

    // RADIANCES---Don'way for indexing
    start4[0] = 0;
    start4[1] = sline;
    start4[2] = 0;
    count2[0] = num_blue_bands;
    count2[1] = 1;  // 1 line at a time
    count2[2] = num_pixels;

    status = nc_get_vara_float(observationGrp, Lt_blueId, start4, count2,
                               l1cstr->Lt_blue[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    count2[0] = num_red_bands;
    status = nc_get_vara_float(observationGrp, Lt_redId, start4, count2, l1cstr->Lt_red[0]);
    check_err(status, __LINE__, __FILE__);

    count2[0] = num_SWIR_bands;
    status = nc_get_vara_float(observationGrp, Lt_SWIRId, start4, count2, l1cstr->Lt_SWIR[0]);
    check_err(status, __LINE__, __FILE__);

    double scantime = ymds2unix((int16_t)scan_time_year, (int16_t)scan_time_month, (int16_t)scan_time_day,
                                scan_time[sline]);
    int16_t syear, sday;
    double secs;
    unix2yds(scantime, &syear, &sday, &secs);

    int32_t yr = syear;
    int32_t dy = sday;
    int32_t msec = (int32_t)(secs * 1000.0);
    double esdist = esdist_(&yr, &dy, &msec);

    float fsol = pow(1.0 / esdist, 2);

    for (size_t ip = 0; ip < num_pixels; ip++) {
        size_t band;
        int ib = 0;
        int ipb = ip * tot_num_bands;

        // load up the blue bands skip the last two
        if(num_blue_bands<expected_num_blue_bands){
            expected_num_blue_bands=num_blue_bands;
        }
        else{
            expected_num_blue_bands=expected_num_blue_bands-2;
        }

        for (band = 0; band < expected_num_blue_bands; band++) {
            if (Lt_blue[band][ip] == Lt_blueFillValue) {
                l1cstr->Lt_tot[ipb] = 0.001;  // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1cstr->Lt_tot[ipb] = Lt_blue[band][ip];
                if (normalizedLt) {
                    l1cstr->Lt_tot[ipb] *= 100;
                } else if (use_rhot) {
                    l1cstr->Lt_tot[ipb] *= l1cstr->Fobar[ib] * fsol * cos(l1cstr->solz[ip] / RADEG) / M_PI;

                } else {
                    l1cstr->Lt_tot[ipb] /= 10.;  // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the red bands skipping the first two and the last 4
        if(num_red_bands<expected_num_red_bands){
            expected_num_red_bands=num_red_bands;
        }
        else{
            expected_num_red_bands=expected_num_red_bands-4;
        }

        for (band = 2; band < expected_num_red_bands; band++) {
            if (Lt_red[band][ip] == Lt_redFillValue) {
                l1cstr->Lt_tot[ipb] = 0.001;  // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1cstr->Lt_tot[ipb] = Lt_red[band][ip];
                if (normalizedLt) {
                    l1cstr->Lt_tot[ipb] *= 100;
                } else if (use_rhot) {
                    l1cstr->Lt_tot[ipb] *= l1cstr->Fobar[ib] * fsol * cos(l1cstr->solz[ip] / RADEG) / M_PI;
                } else {
                    l1cstr->Lt_tot[ipb] /= 10.;  // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }

        // load up the SWIR bands, skip band 3 and 6, hi/low gain wavelengths

        if(num_SWIR_bands<expected_num_SWIR_bands){
            expected_num_SWIR_bands=num_SWIR_bands;
        }
        for (band = 0; band < expected_num_SWIR_bands; band++) {
            if (band == 3 || band == 6)
                continue;
            if (Lt_SWIR[band][ip] == Lt_SWIRFillValue) {
                l1cstr->Lt_tot[ipb] = 0.001;  // should be BAD_FLT, but that makes atmocor fail
            } else {
                l1cstr->Lt_tot[ipb] = Lt_SWIR[band][ip];
                if (normalizedLt) {
                    l1cstr->Lt_tot[ipb] *= 100;
                } else if (use_rhot) {
                    l1cstr->Lt_tot[ipb] *= l1cstr->Fobar[ib] * fsol * cos(l1cstr->solz[ip] / RADEG) / M_PI;
                } else {
                    l1cstr->Lt_tot[ipb] /= 10.;  // the input is in W/m2 ...
                }
            }
            ib++;
            ipb++;
        }
    }

    return 0;
}

int32_t l1c_str::readl1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, int32_t sline) {
    size_t oneline;
    size_t start3[] = {0, 0, 0};
    size_t count2[] = {1, 1, 1};
    string str;

    str = l1cfile->l1b_name;

    printf("Reading HARP2 L1B file\n");
    cout << str << endl;

    // assigning mem for Lt of different bands **********************

    // GEOLOCATION
    l1cstr->iscan = sline;

    // RADIANCES
    start3[0] = 0;
    start3[1] = sline;
    start3[2] = 0;
    count2[0] = nviews;
    count2[1] = 1;  // 1 line at a time
    count2[2] = num_pixels;

    status = nc_get_vara_float(blueGrp, senaId, start3, count2,
                               l1cstr->senazpix_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(blueGrp, latId, start3, count2,
                               l1cstr->latpix_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(blueGrp, lonId, start3, count2,
                               l1cstr->lonpix_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    oneline = sline;

    if (oneline < num_scans - 1)
        start3[1] = sline + 1;
    else
        start3[1] = sline;

    status = nc_get_vara_float(blueGrp, latId, start3, count2,
                               l1cstr->latpix2_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(blueGrp, lonId, start3, count2,
                               l1cstr->lonpix2_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    status = nc_get_vara_float(blueGrp, IntId, start3, count2,
                               l1cstr->I[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    return 0;
}

int32_t l1c_str::readl1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file,
                                  int32_t sline) {
    size_t oneline;
    size_t start3[] = {0, 0, 0};
    size_t count2[] = {1, 1, 1};
    string str;

    str = l1cfile->l1b_name;

    printf("Reading HARP2 L1B file\n");
    cout << str << endl;

    // assigning mem for Lt of different bands **********************

    // GEOLOCATION
    l1cstr->iscan = sline;

    // RADIANCES
    start3[0] = 0;
    start3[1] = sline;
    start3[2] = 0;
    count2[0] = nviews;
    count2[1] = 1;  // 1 line at a time
    count2[2] = num_pixels;

    status = nc_get_vara_float(blueGrp, senaId, start3, count2,
                               l1cstr->senazpix_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(blueGrp, latId, start3, count2,
                               l1cstr->latpix_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(blueGrp, lonId, start3, count2,
                               l1cstr->lonpix_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    oneline = sline;

    if (oneline < num_scans - 1)
        start3[1] = sline + 1;
    else
        start3[1] = sline;

    status = nc_get_vara_float(blueGrp, latId, start3, count2,
                               l1cstr->latpix2_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(blueGrp, lonId, start3, count2,
                               l1cstr->lonpix2_3d[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    status = nc_get_vara_float(blueGrp, IntId, start3, count2,
                               l1cstr->I[0]);  // these are 2-D arrays------------ #bands x #pixels
    check_err(status, __LINE__, __FILE__);

    return 0;
}

//--------------------------------------------------------------------------
//__________________________________________________________________________

int32_t l1c_str::closel1b_spex_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    string str;

    str = l1cfile->l1b_name;
    printf("Closing spexone l1b file\n");
    cout << str << endl;
    status = nc_close(l1cfile->sd_id);
    check_err(status, __LINE__, __FILE__);

    if (l1cstr->viewport != nullptr)
        free (l1cstr->viewport);
    if (l1cstr->senazpix != nullptr)
        free (l1cstr->senazpix);
    if (l1cstr->latpix != nullptr)
        free (l1cstr->latpix);
    if (l1cstr->lonpix != nullptr)
        free (l1cstr->lonpix);
    if (l1cstr->latpix2 != nullptr)
        free (l1cstr->latpix2);
    if (l1cstr->lonpix2 != nullptr)
        free (l1cstr->lonpix2);
    if (l1cstr->latpix3 != nullptr)
        free (l1cstr->latpix3);
    if (l1cstr->lonpix3 != nullptr)
        free (l1cstr->lonpix3);

    if (I_lambdas)
        free(I_lambdas);
    if (pol_lambdas)
        free(pol_lambdas);

    if (l1cstr->I)
        free2d_float(l1cstr->I);
    if (l1cstr->I_polsample)
        free2d_float(l1cstr->I_polsample);

    l1cstr->latpix = nullptr;
    l1cstr->latpix2 = nullptr;
    l1cstr->latpix3 = nullptr;
    l1cstr->lonpix = nullptr;
    l1cstr->lonpix2 = nullptr;
    l1cstr->lonpix3 = nullptr;

    l1cstr->I = nullptr;
    l1cstr->I_polsample = nullptr;

    return 0;
}

int32_t l1c_str::closel1b_oci_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    string str;

    str = l1cfile->l1b_name;
    printf("Closing ocis l1b file\n");
    cout << str << endl;
    status = nc_close(l1cfile->sd_id);
    check_err(status, __LINE__, __FILE__);

    // Free memory
    // From openl1b_ocis
    if (l1cstr->latpix != nullptr)
        free (l1cstr->latpix);
    if (l1cstr->lonpix != nullptr)
        free (l1cstr->lonpix);
    if (l1cstr->latpix2 != nullptr)
        free (l1cstr->latpix2);
    if (l1cstr->lonpix2 != nullptr)
        free (l1cstr->lonpix2);
    if (l1cstr->latpix3 != nullptr)
        free (l1cstr->latpix3);
    if (l1cstr->lonpix3 != nullptr)
        free (l1cstr->lonpix3);
    if (l1cstr->timepix != nullptr)
        free (l1cstr->timepix);
    if (l1cstr->tilt != nullptr)
        free (l1cstr->tilt);

    if (tmpShort)
        free(tmpShort);
    if (scan_time)
        free(scan_time);
    if (blue_lambdas)
        free(blue_lambdas);
    if (red_lambdas)
        free(red_lambdas);
    if (SWIR_lambdas)
        free(SWIR_lambdas);

    if (l1cstr->Lt_blue)
        free2d_float(l1cstr->Lt_blue);
    if (l1cstr->Lt_red)
        free2d_float(l1cstr->Lt_red);
    if (l1cstr->Lt_SWIR)
        free2d_float(l1cstr->Lt_SWIR);

    if (l1cstr->Lt_tot != nullptr)
        free (l1cstr->Lt_tot);
    if (l1cstr->Fobar != nullptr)
        free (l1cstr->Fobar);
    if (l1cstr->solz != nullptr)
        free (l1cstr->solz);
    if (l1cstr->senz != nullptr)
        free (l1cstr->senz);

    l1cstr->latpix = nullptr;
    l1cstr->latpix2 = nullptr;
    l1cstr->latpix3 = nullptr;
    l1cstr->lonpix = nullptr;
    l1cstr->lonpix2 = nullptr;
    l1cstr->lonpix3 = nullptr;
    l1cstr->timepix = nullptr;
    l1cstr->tilt = nullptr;

    l1cstr->Lt_blue = nullptr;
    l1cstr->Lt_red = nullptr;
    l1cstr->Lt_SWIR = nullptr;
    l1cstr->Lt_tot = nullptr;
    l1cstr->Fobar = nullptr;

    return 0;
}

int32_t l1c_str::closel1b_ocis_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    string str;

    str = l1cfile->l1b_name;
    printf("Closing ocis l1b file\n");
    cout << str << endl;
    status = nc_close(l1cfile->sd_id);
    check_err(status, __LINE__, __FILE__);

    // Free memory
    // From openl1b_ocis
    if (l1cstr->latpix != nullptr)
        free (l1cstr->latpix);
    if (l1cstr->lonpix != nullptr)
        free (l1cstr->lonpix);
    if (l1cstr->latpix2 != nullptr)
        free (l1cstr->latpix2);
    if (l1cstr->lonpix2 != nullptr)
        free (l1cstr->lonpix2);
    if (l1cstr->latpix3 != nullptr)
        free (l1cstr->latpix3);
    if (l1cstr->lonpix3 != nullptr)
        free (l1cstr->lonpix3);
    if (l1cstr->timepix != nullptr)
        free (l1cstr->timepix);

    if (tmpShort)
        free(tmpShort);
    if (scan_time)
        free(scan_time);
    if (blue_lambdas)
        free(blue_lambdas);
    if (red_lambdas)
        free(red_lambdas);
    if (SWIR_lambdas)
        free(SWIR_lambdas);

    if (l1cstr->Lt_blue)
        free2d_float(l1cstr->Lt_blue);
    if (l1cstr->Lt_red)
        free2d_float(l1cstr->Lt_red);
    if (l1cstr->Lt_SWIR)
        free2d_float(l1cstr->Lt_SWIR);

    if (l1cstr->Lt_tot != nullptr)
        free (l1cstr->Lt_tot);
    if (l1cstr->Fobar != nullptr)
        free (l1cstr->Fobar);
    if (l1cstr->solz != nullptr)
        free (l1cstr->solz);
    if (l1cstr->senz != nullptr)
        free (l1cstr->senz);

    l1cstr->latpix = nullptr;
    l1cstr->latpix2 = nullptr;
    l1cstr->lonpix = nullptr;
    l1cstr->lonpix2 = nullptr;

    l1cstr->Lt_blue = nullptr;
    l1cstr->Lt_red = nullptr;
    l1cstr->Lt_SWIR = nullptr;
    l1cstr->Lt_tot = nullptr;
    l1cstr->Fobar = nullptr;

    return 0;
}

int32_t l1c_str::closel1b_harp2_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile) {
    string str;

    str = l1cfile->l1b_name;
    printf("Closing harp2 l1b file\n");
    cout << str << endl;
    status = nc_close(l1cfile->sd_id);
    check_err(status, __LINE__, __FILE__);

    if (l1cstr->senazpix_3d)
        free2d_float(l1cstr->senazpix_3d);
    if (l1cstr->latpix_3d)
        free2d_float(l1cstr->latpix_3d);
    if (l1cstr->lonpix_3d)
        free2d_float(l1cstr->lonpix_3d);
    if (l1cstr->latpix2_3d)
        free2d_float(l1cstr->latpix2_3d);
    if (l1cstr->lonpix2_3d)
        free2d_float(l1cstr->lonpix2_3d);

    if (I_lambdas)
        free(I_lambdas);
    if (pol_lambdas)
        free(pol_lambdas);

    if (l1cstr->I)
        free2d_float(l1cstr->I);
    if (l1cstr->I_polsample)
        free2d_float(l1cstr->I_polsample);

    l1cstr->latpix = nullptr;
    l1cstr->latpix2 = nullptr;
    l1cstr->lonpix = nullptr;
    l1cstr->lonpix2 = nullptr;

    l1cstr->I = nullptr;
    l1cstr->I_polsample = nullptr;

    return 0;
}

int32_t l1c_str::closel1b_misr_l1c(l1c_str *l1cstr, l1c_filehandle *l1cfile, filehandle *l1file) {
    string str;

    str = l1cfile->l1b_name;
    printf("Closing harp2 l1b file\n");
    cout << str << endl;
    status = nc_close(l1cfile->sd_id);
    check_err(status, __LINE__, __FILE__);

    if (l1cstr->senazpix_3d)
        free2d_float(l1cstr->senazpix_3d);
    if (l1cstr->latpix_3d)
        free2d_float(l1cstr->latpix_3d);
    if (l1cstr->lonpix_3d)
        free2d_float(l1cstr->lonpix_3d);
    if (l1cstr->latpix2_3d)
        free2d_float(l1cstr->latpix2_3d);
    if (l1cstr->lonpix2_3d)
        free2d_float(l1cstr->lonpix2_3d);

    if (I_lambdas)
        free(I_lambdas);
    if (pol_lambdas)
        free(pol_lambdas);

    if (l1cstr->I)
        free2d_float(l1cstr->I);
    if (l1cstr->I_polsample)
        free2d_float(l1cstr->I_polsample);

    l1cstr->latpix = nullptr;
    l1cstr->latpix2 = nullptr;
    l1cstr->lonpix = nullptr;
    l1cstr->lonpix2 = nullptr;

    l1cstr->I = nullptr;
    l1cstr->I_polsample = nullptr;

    return 0;
}

}  // namespace l1c
