//  bin_L1C.cpp
//
//
//  Created by Martin Montes on 1/30/2023
#include "l1c_bin.h"
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

// static float *Fobar; // reflectance to radiance conversion factors
// static int extract_pixel_start = 0;

// static short *tmpShort;

// whole file stuff
// static size_t expected_num_blue_bands = 60; //120;
// static size_t expected_num_red_bands = 60; //120;
// static size_t expected_num_SWIR_bands = 9;

// static size_t nviews;
// static size_t nbands;
// static size_t num_scans, num_pixels;

// static size_t tot_num_bands = 123; //238;

// scan line attributes

// static int32_t scan_time_hour,scan_time_min;
// static double scan_time_secs;

// static uint_8 *scan_quality;

// geolocation data

// Observation data

namespace l1c {

bin_L1C::bin_L1C() {
    // global attributes
    num_gridlines = -1;
    nbinx = 1;
    nrec_2D = nullptr;  // row/col
    nrec_3D = nullptr;  // row/col/view

    // OCIS binned
    diff_h2 = nullptr;        // row/col
    sca_3D = nullptr;         // row/col/view
    QC_bitwise_4D = nullptr;  // row/col/view/bands
    QC_3D = nullptr;          // row/col/view
    I_4D = nullptr;           // row/col/view/bands
    I_noise_4D = nullptr;

    // OCIS line by line
    obs_per_view = nullptr;
    QC_bitwise = nullptr;
    QC = nullptr;
    I = nullptr;
    I_noise = nullptr;

    l1file = nullptr;

    // ancill info?
}

bin_L1C::~bin_L1C() {
}

int32_t bin_L1C::open_binvars(bin_L1C *binstr, l1c_filehandle *l1cfile) {
    return 0;
}

int32_t bin_L1C::read_binvars(bin_L1C *binstr, l1c_filehandle *l1cfile, int32_t recnum) {
    return 0;
}

int32_t bin_L1C::close_binvars_l1c(bin_L1C *binstr, l1c_filehandle *l1cfile) {
    return 0;
}

}  // namespace l1c
