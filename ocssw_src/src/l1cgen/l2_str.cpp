//
//  l2_str.cpp
//
//
//  Created by Martin Montes on 2/11/2022
#include "l2_str.h"
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

#include "allocate4d.h"
#include <allocate3d.h>
#include <allocate2d.h>

#include <boost/foreach.hpp>
#include <boost/algorithm/string.hpp>

static short *tmpShort;

// whole file stuff
static size_t expected_number_of_bands = 239;  // OCIS
static size_t nviews;
static size_t num_scans, num_pixels, number_of_bands;
static int ncid_L1B;

// scan line attributes
static double *scan_time;  // seconds of day
static int32_t scan_time_year, scan_time_month, scan_time_day;

// static uint_8 *scan_quality;

// geolocation data
static int geolocationGrp;
static int tiltId, lonId, latId, prodims[345];  // #345 of potential L2 products is from product.h
static float latFillValue = -999.0;
static float lonFillValue = -999.0;

static float tiltmin = 0.0, tiltmax = 0.0;
static float tiltFillValue = -999.0;

// Observation data
static int observationGrp;

// Navigation data
static int status, dimid;

using namespace std;

namespace l1c {

l2_str::l2_str()
    : att_ang{-999.0, -999.0, -999.0}, orb_pos{-999.0, -999.0, -999.0}, orb_vel{-999.0, -999.0, -999.0} {
    // init constructor
    // global attributes
    npix = -1;
    iscan = 0;
    nscan = -1;
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

    l2prod = nullptr;
    nl2prod = -1;
    slopeprod = nullptr;
    offsetprod = nullptr;
    tilt = nullptr;

    // sensor/sun geometry
    senz = nullptr;
    sena = nullptr;
    solz = nullptr;
    sola = nullptr;
    delphi = nullptr;
    scattang = nullptr;

    // geolocation--
    senazpix = nullptr;
    latpix = nullptr;
    lonpix = nullptr;
    latpix2 = nullptr;
    lonpix2 = nullptr;

    l1cfile = nullptr;

    // ancill info?
}

l2_str::~l2_str() {
}

int32_t l2_str::openl2_ocis_l1c(L1C_input *l1cinput, l2_str *l2str, l1c_filehandle *l1cfile,
                                int16_t *file_id) {
    std::string str;
    const char *ptstr;

    string delim1 = ":, ";  // product delimiters
    string l2prod_str = l1cinput->l2prod;
    boost::trim_if(l2prod_str, boost::is_any_of(delim1));
    vector<string> prodparam;
    boost::algorithm::split(prodparam, l2prod_str, boost::is_any_of(delim1));
    if(l1cinput->verbose) cout << "number of L2 products to be processed...................#:..." << prodparam.size() << endl;
    for (size_t iprod = 0; iprod < prodparam.size(); iprod++) {
       if(l1cinput->verbose) cout << "selected L2 ----------- prodparam...." << prodparam[iprod] << endl;
    }

    // nl2prod member of l2_str
    nl2prod = prodparam.size();  // number of selected l2 products

    str = l1cfile->l1b_name;
    ptstr = str.c_str();
    printf("Opening OCIS L2 file\n");
    cout << str << endl;

    status = nc_open(ptstr, NC_NOWRITE, &ncid_L1B);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- %s line %d: nc_open(%s) failed.\n", __FILE__, __LINE__, ptstr);
        return (1);
    }

    // num_scans
    status = nc_inq_dimid(ncid_L1B, "number_of_lines", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading number_of_scans.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1B, dimid, &num_scans);
    l2str->nscan = num_scans;

    // num_pixels
    status = nc_inq_dimid(ncid_L1B, "pixels_per_line", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading num_pixels.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1B, dimid, &num_pixels);
    l2str->npix = num_pixels;

    // number of bands
    status = nc_inq_dimid(ncid_L1B, "number_of_bands", &dimid);
    if (status != NC_NOERR) {
        fprintf(stderr, "-E- Error reading number of_bands.\n");
        exit(EXIT_FAILURE);
    }
    nc_inq_dimlen(ncid_L1B, dimid, &number_of_bands);
    if (number_of_bands < expected_number_of_bands) {
        fprintf(stderr, "-E- Not enough  bands, expecting %d, found %d.\n", (int)expected_number_of_bands,
                (int)number_of_bands);
        exit(EXIT_FAILURE);
    }

    //    if (want_verbose) {
    printf("OCI L2 Npix  :%d Nlines:%d\n", (int)num_pixels, (int)num_scans);
    //    } // want_verbose

    // allocate all of the data
    tmpShort = (short *)calloc(num_pixels, sizeof(short));
    scan_time = (double *)calloc(num_scans, sizeof(double));
    l2str->tilt = (float *)calloc(num_scans, sizeof(float));
    l2str->latpix = (float *)calloc(num_pixels, sizeof(float));
    l2str->lonpix = (float *)calloc(num_pixels, sizeof(float));
    l2str->l2prod = allocate2d_float(nl2prod, num_pixels);

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
    if (status == NC_NOERR) {
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
        status = nc_inq_varid(groupid, "msec", &varid);
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

    // Setup geofile pointers
    status = nc_inq_grp_ncid(ncid_L1B, "navigation_data", &geolocationGrp);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(geolocationGrp, "longitude", &lonId);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_var_fill(geolocationGrp, lonId, NULL, &lonFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(geolocationGrp, "latitude", &latId);
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, latId, NULL, &latFillValue);
    check_err(status, __LINE__, __FILE__);

    status = nc_inq_varid(geolocationGrp, "tilt", &tiltId);  // number of lines
    check_err(status, __LINE__, __FILE__);
    status = nc_inq_var_fill(geolocationGrp, tiltId, NULL, &tiltFillValue);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, tiltId, "valid_min", &tiltmin);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_att_float(geolocationGrp, tiltId, "valid_max", &tiltmax);
    check_err(status, __LINE__, __FILE__);

    // get IDs for the observations
    status = nc_inq_grp_ncid(ncid_L1B, "geophysical_data", &observationGrp);
    check_err(status, __LINE__, __FILE__);

    for (size_t iprod = 0; iprod < nl2prod; iprod++) {
        if(l1cinput->verbose) cout << "getting sds id for product.." << prodparam[iprod].c_str() << endl;
        status = nc_inq_varid(observationGrp, prodparam[iprod].c_str(), &prodims[iprod]);
        check_err(status, __LINE__, __FILE__);
    }

    // slope offset of products
    l2str->slopeprod = (float *)calloc(nl2prod, sizeof(float));
    l2str->offsetprod = (float *)calloc(nl2prod, sizeof(float));

    string ATT_NAME1 = "scale_factor", ATT_NAME2 = "add_offset";
    for (size_t iprod = 0; iprod < nl2prod; iprod++) {
        if (nc_get_att_float(observationGrp, prodims[iprod], ATT_NAME1.c_str(), &l2str->slopeprod[iprod]))
            check_err(status, __LINE__, __FILE__);
        if (nc_get_att_float(observationGrp, prodims[iprod], ATT_NAME2.c_str(), &l2str->offsetprod[iprod]))
            check_err(status, __LINE__, __FILE__);
    }

    // tilt
    status = nc_get_var_float(geolocationGrp, tiltId, l2str->tilt);
    check_err(status, __LINE__, __FILE__);

    // l1cfile class---
    l1cfile->sd_id = ncid_L1B;
    l1cfile->nbands = number_of_bands;
    nviews = 2;
    l1cfile->n_views = nviews;

    cout << "number of bands =" << l1cfile->nbands << endl;

    l1cfile->npix = num_pixels;
    l1cfile->nadpix = (num_pixels - 1) / 2;  // nadir pixel index
    l1cfile->nscan = num_scans;
    l1cfile->ndets = 1;
    l1cfile->terrain_corrected = 1;  // presumed.
    l1cfile->orbit_number = orbit_number;

    return 0;
}

int32_t l2_str::readl2_ocis_l1c(l2_str *l2str, l1c_filehandle *l1cfile, int16_t *file_id, int32_t sline) {
    size_t start[] = {0, 0};
    size_t count[] = {1, 1};

    string str;
    str = l1cfile->l1b_name;
    printf("Reading OCIS L2 file\n");
    cout << str << endl;

    // assigning mem for Lt of different bands **********************

    // GEOLOCATION
    l2str->iscan = sline;
    start[0] = sline;
    start[1] = 0;
    count[0] = 1;
    count[1] = num_pixels;  // 1 line at a time

    status = nc_get_vara_float(geolocationGrp, latId, start, count, l2str->latpix);
    check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(geolocationGrp, lonId, start, count, l2str->lonpix);
    check_err(status, __LINE__, __FILE__);

    // L2 products
    for (size_t iprod = 0; iprod < nl2prod; iprod++) {
        status = nc_get_vara_float(observationGrp, prodims[iprod], start, count, l2str->l2prod[iprod]);
        check_err(status, __LINE__, __FILE__);
    }

    return 0;
}

int32_t l2_str::closel2_ocis_l1c(l2_str *l2str, l1c_filehandle *l1cfile) {
    string str;
    str = l1cfile->l1b_name;

    printf("Closing ocis L2 file\n");
    cout << str << endl;
    status = nc_close(l1cfile->sd_id);
    check_err(status, __LINE__, __FILE__);
    // Free memory

    if (l2str->latpix != nullptr)
        free (l2str->latpix);
    if (l2str->lonpix != nullptr)
        free (l2str->lonpix);
    if (l2str->slopeprod != nullptr)
        free (l2str->slopeprod);
    if (l2str->offsetprod != nullptr)
        free (l2str->offsetprod);
    if (l2str->tilt != nullptr)
        free (l2str->tilt);
    if (l2str->l2prod != nullptr)
        free (l2str->l2prod);

    if (tmpShort)
        free(tmpShort);
    if (scan_time)
        free(scan_time);

    l2str->latpix = nullptr;
    l2str->slopeprod = nullptr;
    l2str->lonpix = nullptr;
    l2str->offsetprod = nullptr;
    l2str->tilt = nullptr;
    return 0;
}

}  // namespace l1c
