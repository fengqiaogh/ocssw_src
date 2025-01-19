/*
 * File:   l1c.h
 * Author: mmontes
 *
 * Created on November 4, 2020, 8:45 AM
 * //  last version 8/15/2022
 */

#ifndef L1C_H
#define L1C_H

#include <stdio.h>
#include <string>
#include <vector>
#include <filetype.h>
#include "l1c_filehandle.h"
#include "l1c_str.h"
#include "l2_str.h"
#include "l1c_input.h"
#include "hawkeye_methods.h"
#include <boost/assign/list_of.hpp>  // for 'list_of()'
#include <boost/assert.hpp>
#include <boost/algorithm/string.hpp>
#include <boost/geometry.hpp>
#include <boost/geometry/geometries/point_xy.hpp>
#include <boost/geometry/geometries/polygon.hpp>
#include <boost/geometry/geometries/box.hpp>
#include <boost/foreach.hpp>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <netcdf>

#define READ 0
#define WRITE 1

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> Point_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::box<Point_t> Box_t;
typedef double orb_array2[3];

namespace l1c {

class L1C {
   protected:
   public:
    // methods--------------------------------------
    L1C();
    virtual ~L1C();

    virtual int32_t load_l1c_filehandle4(l1c_filehandle *l1cfile, L1C_input *l1cinput);

    virtual int32_t ect_swt(l1c_filehandle *l1cfile, int ix1, int ix2, double *tswt_tot, double *latswt_tot,
                            double *lonswt_tot, double *ovel_tot, double *gvel_tot, double *tswt,
                            double *latswt, double *lonswt, float *tcross, float *loncross, double *ovel,
                            double *gvel);
    virtual int32_t create_time_swt(int num_gridlines, double tfile_ini_sec, double *tmgvf,
                                    double tswt_ini_sec, double tswt_end_sec, std::string *tswt_ini,
                                    std::string *tswt_ini_file, std::string *tswt_mid,std::string *tswt_end);
    virtual int32_t swtime_swt2(int swt, L1C_input *l1cinput, l1c_filehandle *l1cfile, int32_t norbs,
                                double *tswt, double tcross, double mgv, double *tmgv, double *orb_time_tot,size_t norbs_tot);
    virtual int32_t swtime_swt2_segment(int swt, L1C_input *l1cinput, l1c_filehandle *l1cfile, int32_t norbs,
                                        double *tswt, double tcross, double mgv, double *tmgv,double *orb_time_tot,size_t norbs_tot);
    virtual int32_t write_L1C_granule2(int swtd, l1c_filehandle *l1cfile, L1C_input *l1cinput, double *tmgv,
                                       float **lat_gd, float **lon_gd, float **alt_gd,double *orb_time_tot);
    virtual int32_t open_l1atol1c3(L1C_input *l1cinput, l1c_filehandle *l1cfile);
    virtual int search_l1cgen(L1C_input *l1cinput, l1c_str *l1cstr, l1c_filehandle *l1cfile, short **gdindex);
    virtual int32_t create_SOCEA2(int swtd, L1C_input *l1cinput, l1c_filehandle *l1cfile, float **lat_gd,
                                  float **lon_gd, float **altitude, double *tswt);
    virtual int32_t openL1Cgrid3(l1c_str *l1cstr, l1c_filehandle *l1cfile, L1C_input *l1cinput);
    virtual int32_t l1_cloud_correct(L1C_input *l1cinput, l1c_filehandle *l1cfile);
    virtual int32_t add_proc_group_l1c(L1C_input *l1cinput,l1c_filehandle *l1cfile,const char* filename);

    //---------------------------------------------
    // global attributes----
    std::string l1b_name;  // just the input l1b filenamed
    size_t sd_id;

    // more L1C input vars
    size_t l1c_pflag;                       // L1C processing
    std::vector<std::string> cust_l1cprod;  // list of L1c products to be included, flexible approach

    file_type format;  // file type, netcdf4 or hdf4/5
    // size_t mode;//l1c processing mode ---
    size_t length;    // data block
    size_t sensorID;  // original int32 changed to size_t type
    size_t subsensorID;
    float res_spat;  // spatial resolution in km
    float res_spec;  // spectral resolution in nm

    // dimensions--
    size_t ndets;
    size_t nscan;
    size_t n_views;     // sensor views
    size_t npols;       // polarization states
    size_t nbands;      // number of total bands
    size_t nband_blue;  // this includes uv + visible bands
    size_t nband_red;   // this inc
    size_t nband_swir;
    size_t npix;  // num

    // sensor characteritics
    float *views;    // indexes are not defined
    size_t pols[3];  //[1 0 0] non-pol, [1 1 0]: non+cross, [1 1 1]: non+cross+parall
    float *bbands;   // array with bands wavelengths
    float *rbands;
    float *swirbands;

    // navigation attributes
    // fred attr here
    size_t orbit_number;
    size_t orb_dir;        // asc or desc orbit
    float orbit_node_lon;  // long at which is crossing the equator asc or desc

    // geolocation attributes
    size_t terrain_corrected;
    size_t cloud_corrected;
    float *cloud_height;  // should be in geolocation group nc file

    // calibration attributes
    float *Fobar;

    // projection attributes-
    // std::string proj_type;
    size_t projection;
    float grid_resolution;  // grid resolution in km
    // these params go to proj class

    // multi  attributes (view, pol, bands)
    float *view_agg;  // views to be aggregated for later products such as vsfs etc
    float *pol_agg;   // polarization states to be aggregated for post-processing products, linear
                     // depolarization ratio etc
    float *band_agg;     // specific bands for future merged products
    bool overlap_vflag;  // tells if we want merged views
    bool overlap_pflag;  // tells if we want merged polarizations
    bool overlap_bflag;  // tells if we want merged spectral bands
    // uncertainty params for merged l1c products
    size_t unc_meth;    // 0: no error calculation, 1: propagation, 2: Monte Carlo
    float unc_thres_v;  // uncertainity threshold of angular merged products as %
    float unc_thres_p;  // same but for polarization
    float unc_thres_b;  // same but for spectral bands
    // ancill info, requested products?
    // float *cloud_h;//this is cloud_height
};

// protos---
int32_t load_l1c_filehandle4(l1c_filehandle *l1cfile, L1C_input *l1cinput);
int32_t ect_swt(l1c_filehandle *l1cfile, int ix1, int ix2, double *tswt_tot, double *latswt_tot,
                double *lonswt_tot, double *ovel_tot, double *gvel_tot, double *tswt, double *latswt,
                double *lonswt, float *tcross, float *loncross, double *ovel, double *gvel);
int32_t create_time_swt(int num_gridlines, double tfile_ini_sec, double *tmgvf, double tswt_ini_sec,
                        double tswt_end_sec, std::string *tswt_ini, std::string *tswt_ini_file,std::string *tswt_mid,
                        std::string *tswt_end);
int32_t swtime_swt2(int swt, L1C_input *l1cinput, l1c_filehandle *l1cfile, int32_t norbs, double *tswt,
                    double tcross, double mgv, double *tmgv,double *orb_time_tot,size_t norbs_tot);
int32_t swtime_swt2_segment(int swt, L1C_input *l1cinput, l1c_filehandle *l1cfile, int32_t norbs,
                            double *tswt, double tcross, double mgv, double *tmgv,double *orb_time_tot,size_t norbs_tot);
int32_t write_L1C_granule2(int swtd, l1c_filehandle *l1cfile, L1C_input *l1cinput, double *tmgv,
                           float **lat_gd, float **lon_gd, float **alt_gd,double *orb_time_tot);
int32_t open_l1atol1c3(L1C_input *l1cinput, l1c_filehandle *l1cfile);
int search_l1cgen(L1C_input *l1cinput, l1c_str *l1cstr, l1c_filehandle *l1cfile, short **gdindex);
int32_t create_SOCEA2(int swtd, L1C_input *l1cinput, l1c_filehandle *l1cfile, float **lat_gd, float **lon_gd,
                      float **altitude, double *tswt);
int32_t openL1Cgrid3(l1c_str *l1cstr, l1c_filehandle *l1cfile, L1C_input *l1cinput);
int32_t add_proc_group_l1c(L1C_input *l1cinput,l1c_filehandle *l1cfile,const char* filename);

}  // namespace l1c
#endif
