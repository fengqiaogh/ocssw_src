
//  l1c_filehandle.cpp
//  Created by Martin Montes on 8/15/22

#include "l1c_filehandle.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <boost/assign/list_of.hpp>  // for 'list_of()'
#include <boost/assert.hpp>
#include <list>
#include <stack>
#include <vector>
#include <iostream>

using namespace std;
using namespace boost::assign;

namespace l1c {

l1c_filehandle::l1c_filehandle() {
    // global--
    // L1C produts-
    // std::vector<std::string> cust_l1cprod = list_of("pc")("vsf")("dpr");//3 options selected for  for l1c
    // products, principal components, volume scattering function and degree linear polarization
    // std::vector<std::string> cust_l1cprod = {"pc","vsf","dpr"};
    cust_l1cprod.push_back("pc");
    gridname = "";
    azeast_name = "";
    swath_scans = 1;
    swath_num = 1;
    for (int j = 0; j < 10; j++) {
        selgran[j] = -1;
    }
    verbose=0;
    l1c_pflag = 0;
    l1b_name = "";
    version = "";
    progname = "";
    sd_id = 0;
    format = FT_INVALID;
    // mode=0;
    length = 0;
    sensorID = -1;
    subsensorID = -1;
    res_spat = -999.0;
    ;
    res_spec = -999.0;
    lat_gd = nullptr;
    lon_gd = nullptr;
    alt_gd = nullptr;
    lat_asort = nullptr;
    index_xy = nullptr;

    // time-space limits of image
    syear = 0;
    sday = 0;
    minlat_img = -999.0;
    ;
    maxlat_img = -999.0;
    ;
    minlon_img = -999.0;
    ;
    maxlon_img = -999.0;
    ;

    // dimensions--
    nframes = -1;  // HARP sensor
    ndets = 1;
    nscan = 0;
    n_views = 1;     // sensor views
    npols = 1;       // polarization states
    nbands = 0;      // number of total bands
    nband_blue = 0;  // this includes uv + visible bands
    nband_red = 0;   // this includes nir + visible bands
    nband_swir = 0;
    // spex------------
    nband_view = 0;
    npol_band_view = 0;
    npix = 0;

    // sensor characteritics
    views = nullptr;  // indexes are not defined
    for (int i = 0; i < 3; i++) {
        pols[i] = 0;  // 0: non-pol, 1: cross, 2: parall
    }

    bbands = nullptr;  // array with bands wavelengths
    rbands = nullptr;
    swirbands = nullptr;

    // navigation attributes
    orbit_number = 0;
    orb_dir = -1;             // 0 asc 1 des
    orbit_node_lon = -999.0;  // this is node_crossing_time
    eqt = -999.0;
    tswt_ini = "";  // metadata string for initial time
    tswt_end = "";  // metadata string for final time
    tswt_mid="";
    tswt_ini_file = "";
    tunix_start = -999.;  // start utc time for the swath,
    tg_s = -999.;
    tg_e = -999.;
    numgran = -1;
    tfile_ini_sec = -999.;
    tfile_end_sec = -999.;
    instrument = "";
    start_dir = "";
    end_dir = "";
    binstr = "";
    date_created = "";
    norb_rec=-1;//number of orbital records total for swath
    // telemetry stauff----L1A--HKT
    fileix = 0;
    gransize = 5;  // granule size in minutes
    gd_per_gran = 400;
    gran_ini_timeflag = 1;  // initial time for granules produced HKT->L1C

    // geolocation attributes
    terrain_corrected = 0;
    cloud_corrected = 0;
    cloud_height = nullptr;

    // calibration attributes
    Fobar = nullptr;

    // projection attribute
    mean_az_east = -999.0;
    NY1 = -1;
    NY2 = -1;
    num_gridlines = -1;
    binyb = -1;
    binya = -1;
    mot = -999.;
    nbinx = -1;
    sensor = -1;
    latnorth = 90.0;
    latsouth = -90.0;
    lonwest = -180.;
    loneast = +180.;
    nswath = -1;
    ndswaths = -1;
    nswt_files = -1;
    nadpix = -1;
    gres = -999.0;

    proj_type = 0;  // 0 'socea" and 1: socea-2
    sort_type = 0;
    lat0 = 0.0;

    // multi  attributes (view, pol, bands)
    view_agg = nullptr;  // views to be aggregated for later products such as vsfs etc
    pol_agg = nullptr;   // polarization states to be aggregated for post-processing products, linear
                        // depolarization ratio etc
    band_agg = nullptr;  // specific bands for future merged products
    overlap_vflag = 0;   // tells if we want merged views
    overlap_pflag = 0;   // tells if we want merged polarizations
    overlap_bflag = 0;   // tells if we want merged spectral bands
    // uncertainty params l1c merged products
    unc_meth = -1;
    unc_thres_v = -999.0;  // uncertainity threshold of angular merged products as %
    unc_thres_p = -999.0;  // same but for polarization
    unc_thres_b = -999.0;  // same but for spectral bands
}

l1c_filehandle::~l1c_filehandle() {
}

void l1c_filehandle::printProducts() {
    std::vector<std::string> cust_l1cprod = {"pc", "vsf", "dpr"};
    // cust_l1cprod = {"pc", "vsf", "dpr"};
    for (const auto& item : cust_l1cprod) {
        std::cout << item << std::endl;
    }
}

}  // namespace l1c
