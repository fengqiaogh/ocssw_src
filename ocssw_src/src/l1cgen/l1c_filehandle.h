/*
 * File:   l1c_filehandle.h
 * Author: mmontes
 *
---------------------------------------------------------
l1cgen firstly defined as mainL1C.cpp
  Created by Martin Montes on 8/15/2022
  last version 3/7/2022
---------------------------------------------------------
 */

#ifndef L1C_FILEHANDLE_H
#define L1C_FILEHANDLE_H

#include <stdio.h>
#include <string>
#include <vector>
#include <filetype.h>

#define READ 0
#define WRITE 1

namespace l1c {

class l1c_filehandle {
   protected:
   public:
    // methods---
    l1c_filehandle();
    virtual ~l1c_filehandle();
    virtual void printProducts();
    // global attributes----
    std::string version;
    std::string progname;
    std::string l1b_name;  // just the input l1b filenamed
    size_t sd_id;

    // more L1C input vars
    bool verbose;
    size_t swath_num;
    int16_t swath_scans;
    int16_t selgran[10];      // selected granule ids for L1C processing
    const char *gridname;     // name for l1cgrid having 2-d lon,lat gd points
    const char *azeast_name;  // filename of az_east in degrees
    size_t l1c_pflag;  // L1C processing 0:no processing, 1:full orbit, 2: full swath, half-day daytime, 3:
                       // L1C grid generation, 4: 5-minute granule file
    std::vector<std::string> cust_l1cprod;  // list of L1c products to be included, flexible approach
    std::vector<std::string> ifiles;
    file_type format;  // file type, netcdf4 or hdf4/5
    // size_t mode;//read/write file mode
    size_t length;    // data block
    size_t sensorID;  // original int32 changed to size_t type
    size_t subsensorID;
    float res_spat;  // spatial resolution in km
    float res_spec;  // spectral resolution in nm

    // time-space limits of image
    int16_t syear;
    int16_t sday;
    float minlat_img;
    float maxlat_img;
    float minlon_img;
    float maxlon_img;

    // dimensions--
    size_t nframes;
    size_t ndets;
    size_t nscan;
    size_t n_views;     // sensor views
    size_t npols;       // polarization states
    size_t nbands;      // number of total intensity bands
    size_t nband_blue;  // this includes uv + visible bands
    size_t nband_red;   // this includes nir + visible bands
    size_t nband_swir;
    size_t nband_view;
    size_t npol_band_view;
    size_t npix;  // number of pixels per line or across-track --x component

    // sensor characteritics
    float *views;    // indexes are not defined
    size_t pols[3];  //[1 0 0] non-pol, [1 1 0]: non+cross, [1 1 1]: non+cross+parall
    float *bbands;   // array with bands wavelengths
    float *rbands;
    float *swirbands;

    // navigation attributes
    // swath attr here
    size_t orbit_number;
    size_t orb_dir;        // asc or desc orbit 1 or 0 respectively
    float orbit_node_lon;  // long at which is crossing the equator asc or desc
    double eqt;            // equat crossing time
    std::string tswt_ini;  // UTC initial time of swath
    std::string tswt_end;  // UTC ennd time of swath
    std::string tswt_mid;
    std::string tswt_ini_file;
    std::string start_dir;
    std::string end_dir;
    std::string binstr;
    std::string instrument;
    std::string date_created;

    double tunix_start;  // unix time
    double tg_s;
    double tg_e;
    int numgran;
    double tfile_ini_sec;  // unix time seconds since 1970
    double tfile_end_sec;

    // telemetry stauff----L1A--HKT
    int fileix;
    int gransize;           // granule size in minutes
    int32_t gd_per_gran;    // #gridlines per granule time
    int gran_ini_timeflag;  // initial time for L1C granule

    // geolocation attributes
    size_t terrain_corrected;
    size_t cloud_corrected;
    float *cloud_height;  // should be in geolocation group nc file

    // calibration attributes
    float *Fobar;

    // L1C projection attributes-
    // std::string proj_type;
    float mean_az_east;  // swath mean bearing for right side scan
    float **lat_gd;
    float **lon_gd;
    float **alt_gd;
    float **lat_asort;
    short **index_xy;
    int32_t NY1;
    int32_t NY2;
    int32_t num_gridlines;
    int16_t binyb;
    int16_t binya;
    double mot;
    int32_t nbinx;  // grid bins across-track
    int32_t sensor;
    float latnorth;
    float latsouth;
    float lonwest;
    float loneast;
    int16_t nswath;
    int16_t ndswaths;
    int16_t nswt_files;
    int32_t nadpix;
    size_t proj_type;
    size_t sort_type;
    float gres;  // grid resolution in km
    size_t norb_rec;
    // these params go to proj class
    float lat0;

    // multi  attributes (view, pol, bands)
    float *view_agg;  // views to be aggregated for later products such as vsfs etc
    float *pol_agg;   // polarization states to be aggregated for post-processing products, linear
                     // depolarization ratio etc
    float *band_agg;       // specific bands for future merged products
    size_t overlap_vflag;  // tells if we want merged views
    size_t overlap_pflag;  // tells if we want merged polarizations
    size_t overlap_bflag;  // tells if we want merged spectral bands
    // uncertainty params for merged l1c products
    size_t unc_meth;    // 0: no error calculation, 1: propagation, 2: Monte Carlo
    float unc_thres_v;  // uncertainity threshold of angular merged products as %
    float unc_thres_p;  // same but for polarization
    float unc_thres_b;  // same but for spectral bands
    // ancill info, requested products?
    // float *cloud_h;//this is cloud_height
};

// we move this to l1c class --void load_l1c_filehandle(filehandle *file,l1c_filehandle *filel1c);
// void printVector(const std::vector<std::string> &n);
void printProducts();

}  // namespace l1c
#endif
