
//  Created by Martin Montes on 10/26/20.
// latest update 6/28/2023
//****************************************************8

#include <filehandle.h>
#include <l1.h>
#include "l1c.h"
#include "l1c_filehandle.h"
#include "l1c_input.h"
#include "l1c_str.h"
#include "l2_str.h"
#include <string>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <list>
#include <stack>
#include <vector>
#include <iostream>
#include <fstream>
#include <sstream>
#include <netcdf.h>
#include <nc4utils.h>
#include "hawkeye_methods.h"
#include "allocate4d.h"
#include <allocate3d.h>
#include <allocate2d.h>
#include <algorithm>
// #include <bits/stdc++.h>
#include <cstdlib>
#include <ctime>
#include <chrono>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/GeodesicLine.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Rhumb.hpp>
#include <timeutils.h>
#include <genutils.h>
#include <stdint.h>
#include "l2prod.h"
#include "Strassen.h"
#include <libnav.h>

static size_t num_scans, num_pixels;
static size_t num_blue_bands, num_red_bands, num_SWIR_bands;

// Navigation data
static short **senz = nullptr, **sena = nullptr;

// static float *orbp;
// static float *attang;
// static double *timeArr;

// L1C vars
static int16_t nswath;
float  Re = 6378.137;
using namespace std;
using namespace boost::assign;
using namespace std::chrono;
using namespace GeographicLib;

using namespace netCDF;
using namespace netCDF::exceptions;

namespace bg = boost::geometry;
typedef bg::model::point<double, 2, bg::cs::geographic<bg::degree>> Point_t;
typedef bg::model::polygon<Point_t> Polygon_t;
typedef bg::model::box<Point_t> Box_t;

namespace l1c {

L1C::L1C() {
}

L1C::~L1C() {
}

int L1C::add_proc_group_l1c(L1C_input *l1cinput,l1c_filehandle *l1cfile,const char* l1c_filename)
{
//additional attributes

   int ncid_out,grp_pc,grp_ip;
   string str,str2;
   int status = nc_open(l1c_filename, NC_WRITE, &ncid_out);
   if ((status = nc_def_grp(ncid_out, "processing_control", &grp_pc)))  // netcdf-4
            check_err(status, __LINE__, __FILE__);
   str="software_name";
   str2=l1cfile->progname;
   status = nc_put_att_text(grp_pc,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
   str="software_version";
   str2=l1cfile->version;
   str2=str2.substr(0,5);
   status = nc_put_att_text(grp_pc,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());

   if ((status = nc_def_grp(grp_pc, "input_parameters", &grp_ip)))  // netcdf-4
            check_err(status, __LINE__, __FILE__);
   str="ifile";
   str2="";
   int nfiles=l1cfile->ifiles.size();
      for(int f=0;f<nfiles;f++)
      {
      if(f==nfiles-1)str2=str2+l1cfile->ifiles[f];else str2=str2+l1cfile->ifiles[f]+",";
                   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
      }
   str="ofile";
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),strlen(l1c_filename),l1c_filename);
   str="outlist";
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),strlen(l1cinput->outlist),l1cinput->outlist);
   str="l1c_grid";
   str2=string(l1cinput->l1c_grid);
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
   str="l1c_anc";
   str2=string(l1cinput->l1c_anc);
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
   str="demfile";
   str2="$OCDATAROOT/common/gebco_ocssw_v2020.nc";
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());

   str="verbose";
   int value=(int)l1cinput->verbose;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="mode";
   value=l1cinput->l1c_pflag;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="north";
   float value2=l1cinput->north;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="south";
   value2=l1cinput->south;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="east";
   value2=l1cinput->east;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="west";
   value2=l1cinput->west;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="selgran";
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,10,l1cinput->selgran);
   str="selday";
   value=l1cinput->selday;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="selmon";
   value=l1cinput->selmon;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="selyear";
   value=l1cinput->selyear;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="gransize";
   value=l1cinput->gransize;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="grantype";
   value=l1cinput->grantype;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="bintype";
   value=l1cinput->bintype;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="start_timeflag";
   value=l1cinput->start_timeflag;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="swath_num";
   value=l1cinput->swath_num;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="grid_resolution";
   value2=l1cinput->grid_resolution;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="sensor";
   value=l1cinput->sensor;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="l2prod";
   str2=l1cinput->l2prod;
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
   str="ix_l1cprod";
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,3,l1cinput->ix_l1cprod);
   str="start_time";
   str2=l1cinput->start_time;
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
   str="end_time";
   str2=l1cinput->end_time;
   status = nc_put_att_text(grp_ip,NC_GLOBAL,str.c_str(),str2.size(),str2.c_str());
   str="projection";
   value=l1cinput->projection;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="sort_method";
   value=l1cinput->sort_method;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="demcloud_flag";
   value=l1cinput->demcloud_flag;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="terrain_correct";
   value=(int)l1cinput->terrain_correct;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="cloud_height";
   value=l1cinput->cloud_height;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="cloud_correct";
   value=l1cinput->cloud_correct;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="cloud_type";
   value=l1cinput->cloud_type;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="overlap_vflag";
   value=(int)l1cinput->overlap_vflag;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="overlap_pflag";
   value=(int)l1cinput->overlap_pflag;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="overlap_bflag";
   value=(int)l1cinput->overlap_bflag;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="unc_meth";
   value=l1cinput->unc_meth;
   status = nc_put_att_int(grp_ip,NC_GLOBAL,str.c_str(),NC_INT,1,&value);
   str="unc_thres_v";
   value2=l1cinput->unc_thres_v;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="unc_thres_p";
   value2=l1cinput->unc_thres_p;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   str="unc_thres_b";
   value2=l1cinput->unc_thres_b;
   status = nc_put_att_float(grp_ip,NC_GLOBAL,str.c_str(),NC_FLOAT,1,&value2);
   status = nc_close(ncid_out);

return status;
}

// function to calculate cross product of two vectors, 3 components each vector
void cross_product_double(double vector_a[], double vector_b[], double temp[]) {
    temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
    temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];
}

double cross_product_norm_double(double vector_a[], double vector_b[]) {
    double temp[3], nvec;
    temp[0] = vector_a[1] * vector_b[2] - vector_a[2] * vector_b[1];
    temp[1] = -(vector_a[0] * vector_b[2] - vector_a[2] * vector_b[0]);
    temp[2] = vector_a[0] * vector_b[1] - vector_a[1] * vector_b[0];

    nvec = sqrt(temp[0] * temp[0] + temp[1] * temp[1] + temp[2] * temp[2]);
    return nvec;
}

int32_t L1C::l1_cloud_correct(L1C_input* l1cinput, l1c_filehandle* l1cfile) {
    double *ptime = nullptr, *ptime_l1c = nullptr;
    const char* l1c_gran;
    size_t num_ybin, num_xbin;
    float **latpix = nullptr, **lonpix = nullptr, **l1clat = nullptr, **l1clon = nullptr, **l1calt = nullptr;
    float **lat_new = nullptr, **lon_new = nullptr;
    float Rpole = 6356;
    float Xcorr,Ycorr, Zcorr;
    float Rsurf, radius_ratio = Re / Rpole;
    float latcorr;
    size_t gran = 0;
    float dv, Hsat = 676.5;  // sensor height in meters above ellipsoid
    short **senapix = nullptr, **senzpix = nullptr;
    float **cth = nullptr, temp;
    NcFile *nc_l1c, *nc_l1c_new, *nc_l12_new, *nc_l12;
    string name;
    double FillValue2=BAD_FLT;
    float FillValue=BAD_FLT;

    if (l1cinput->cloud_correct == 0) {
        cout << "cloud_correct is 0,thus no paralallax correction is needed!!..use 1 (constant cloud height) "
                "or 2 (variable cloud height)....."
             << endl;
        exit(1);
    } else if (l1cinput->cloud_correct == 1) {
        cout << "l1cinput->cloud_correct is: " << l1cinput->cloud_correct
             << ".LIC->(L1C paralallax corrected)...cloud height in km is constant throughout the L1C grid "
                "but different from zero so parallax correction id needed..."
             << endl;
        if (l1cinput->cloud_height == 0) {
            cout << "ERROR cloud_height is NOT greater than 0..or was not provided as input:"
                 << l1cinput->cloud_height << endl;
            exit(1);
        }

        //--- OPEN L1C FILES ---
        for (unsigned int i = 0; i < l1cinput->files_l1c.size(); i++) {
            gran++;
            string l1c_str =
                l1cinput->files_l1c[i];  // ifile is used here for old grid with no parallax correction, this
                                         // is a list of l1c grid granules or full grid L1C
            l1c_gran = l1c_str.c_str();

            cout << "Opening L1C granule .#.." << i + 1 << "...." << l1c_str
                 << "for parallax-cloud corrections......." << endl;

            try {
                nc_l1c = new NcFile(l1c_gran, NcFile::read);
            } catch (NcException& e) {
                e.what();
                cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: " + l1c_str << endl;
                exit(1);
            }

            NcDim yd = nc_l1c->getDim("bins_along_track");
            num_ybin = yd.getSize();
            NcDim xd = nc_l1c->getDim("bins_across_track");
            num_xbin = xd.getSize();

            ptime_l1c = (double*)calloc(num_ybin, sizeof(double));
            l1clat = allocate2d_float(num_ybin, num_xbin);
            l1clon = allocate2d_float(num_ybin, num_xbin);
            l1calt = allocate2d_float(num_ybin, num_xbin);

            NcGroup ba_grp = nc_l1c->getGroup("bin_attributes");
            NcVar v1 = ba_grp.getVar("nadir_view_time");
            v1.getVar(ptime_l1c);

            NcGroup geo_grp = nc_l1c->getGroup("geolocation_data");
            v1 = geo_grp.getVar("latitude");
            v1.getVar(&l1clat[0][0]);
            v1 = geo_grp.getVar("longitude");
            v1.getVar(&l1clon[0][0]);
            v1 = geo_grp.getVar("altitude");
            v1.getVar(&l1calt[0][0]);

            cout << "num along bins L1C..." << num_ybin << "num_across bins L1C...." << num_xbin
                 << "l1c time_ini.." << ptime_l1c[0] << "l1c time end.." << ptime_l1c[num_ybin - 1] << endl;

            lat_new = allocate2d_float(num_ybin, num_xbin);
            lon_new = allocate2d_float(num_ybin, num_xbin);

            for (size_t i = 0; i < num_ybin; i++) {
                for (size_t j = 0; j < num_xbin; j++) {
                    // Rsurf is Earth radius at pixel positions
                    Rsurf = Re / (sqrt(cos(l1clat[i][j] * M_PI / 180) * cos(l1clat[i][j] * M_PI / 180) +
                                       radius_ratio * radius_ratio * sin(l1clat[i][j] * M_PI / 180) *
                                           sin(l1clat[i][j] * M_PI / 180)));
                    Xcorr = (l1cinput->cloud_height + Rsurf) * cos(l1clat[i][j] * M_PI / 180) *
                            sin(l1clon[i][j] * M_PI / 180);  // lon speherical and geodetic should be equal
                    Ycorr = (l1cinput->cloud_height + Rsurf) * sin(l1clat[i][j] * M_PI / 180);
                    Zcorr = (l1cinput->cloud_height + Rsurf) * cos(l1clat[i][j] * M_PI / 180) *
                            cos(l1clon[i][j] * M_PI / 180);
                    //     convert back to latitude and longitude
                    latcorr = atan(Ycorr / sqrt(Xcorr * Xcorr + Zcorr * Zcorr));
                    lat_new[i][j] = atan(tan(latcorr) / radius_ratio * radius_ratio) * 180.0 / M_PI;
                    lon_new[i][j] = atan2(Xcorr, Zcorr) * 180.0 / M_PI;
                    lon_new[i][j] = atan2(Ycorr, Xcorr) * 180 / M_PI;

                    if(l1cinput->verbose) cout<<"l1cinput->cloud_height.."<<l1cinput->cloud_height<<"latnew.."<<lat_new[i][j]<<"lon_new.."<<lon_new[i][j]<<"l1clat.."<<l1clat[i][j]<<"l1clon.."<<l1clon[i][j]<<endl;
                }
            }

            // write the new L1C gridor grids
            string gran_str = to_string(gran);
            string cth_str = to_string(l1cinput->cloud_height);
            string l1c_gran_new = l1c_str.substr(0, 24) + "_CTH" + cth_str + ".gran" + gran_str + ".nc";
            if(l1cinput->verbose)cout << "corrected L1C grid" << l1c_gran_new << "     granule #" << gran_str << endl;

            const char* l1c_gran = l1c_gran_new.c_str();

            try {
                nc_l1c_new = new NcFile(l1c_gran, NcFile::replace);
            } catch (NcException& e) {
                e.what();
                cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: " + l1c_gran_new << endl;
                exit(1);
            }
            // global attributes--
            char* gridchar = strdup(l1c_gran);

            bin_str binl1c;
            meta_l1c_grid(gridchar, &binl1c, num_ybin, nc_l1c_new);

            nc_l1c_new->putAtt("processing_version", l1cinput->pversion);
            nc_l1c_new->putAtt("history", l1cinput->history);
            string name;
            nc_l1c_new->putAtt("product_name", l1c_gran_new);
            NcGroupAtt i1 = nc_l1c->getAtt("startDirection");  // NcGroupAtt is a global attr!!
            i1.getValues(name);
            nc_l1c_new->putAtt("startDirection", name);
            i1 = nc_l1c->getAtt("endDirection");
            i1.getValues(name);
            nc_l1c_new->putAtt("endDirection", name);
            i1 = nc_l1c->getAtt("time_coverage_start");
            i1.getValues(name);
            nc_l1c_new->putAtt("time_coverage_start", name);
            i1 = nc_l1c->getAtt("time_coverage_end");
            i1.getValues(name);
            nc_l1c_new->putAtt("time_coverage_end", name);
            i1 = nc_l1c->getAtt("date_created");
            i1.getValues(name);
            nc_l1c_new->putAtt("date_created", name);

            ba_grp = nc_l1c_new->getGroup("bin_attributes");
            v1 = ba_grp.getVar("nadir_view_time");
            v1.putVar(&ptime_l1c[0]);
            geo_grp = nc_l1c_new->getGroup("geolocation_data");
            v1 = geo_grp.getVar("latitude");
            v1.putVar(&lat_new[0][0]);
            v1 = geo_grp.getVar("longitude");
            v1.putVar(&lon_new[0][0]);
            v1 = geo_grp.getVar("altitude");
            v1.putVar(&l1calt[0][0]);

            nc_l1c->close();
            nc_l1c_new->close();

            free (lat_new);
            free (lon_new);
            free (ptime_l1c);
            free (l1clat);
            free (l1clon);
            free (l1calt);
        }
    } else if (l1cinput->cloud_correct == 2) {
        if(l1cinput->verbose)cout << "l1cinput->cloud_correct :" << l1cinput->cloud_correct
             << "LIB->(L1B paralallax corrected)----cloud_height is constant between pixels......" << endl;

        if (l1cinput->cloud_height == 0) {
            cout << "ERROR cloud_height is NOT greater than 0..or was not provided as input:"
                 << l1cinput->cloud_height << endl;
            exit(1);
        }

        for (unsigned int i = 0; i < l1cinput->files.size(); i++) {
            gran++;
            string l12_str = l1cinput->files[0];
            const char* l_12 = l12_str.c_str();

            if(l1cinput->verbose)cout << "Opening L1B/L2 granule ..." << l1cinput->files[0]
                 << "....for parallax-cloud corrections......." << endl;

            try {
                nc_l12 = new NcFile(l_12, NcFile::read);
            } catch (NcException& e) {
                e.what();
                cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: " + l12_str << endl;
                exit(1);
            }

            // DIMENSIONS
            NcDim yd = nc_l12->getDim("number_of_scans");
            NcDim xd = nc_l12->getDim("ccd_pixels");
            num_scans = yd.getSize();
            num_pixels = xd.getSize();

            lat_new = allocate2d_float(num_scans, num_pixels);
            lon_new = allocate2d_float(num_scans, num_pixels);

            if(l1cinput->verbose)cout << "number of scans.." << num_scans << "num of pixels.." << num_pixels << endl;

            latpix = allocate2d_float(num_scans, num_pixels);
            lonpix = allocate2d_float(num_scans, num_pixels);
            ptime = (double*)calloc(num_scans, sizeof(double));
            senzpix = allocate2d_short(num_scans, num_pixels);
            senapix = allocate2d_short(num_scans, num_pixels);
            cth = allocate2d_float(num_scans, num_pixels);

            for (unsigned int i = 0; i < num_scans; i++) {
                for (unsigned int j = 0; j < num_pixels; j++) {
                    cth[i][j] = l1cinput->cloud_height;
                }
            }

            NcGroup ba_grp = nc_l12->getGroup("scan_line_attributes");
            NcVar v1 = ba_grp.getVar("time");
            v1.getVar(ptime);
            NcGroup geo_grp = nc_l12->getGroup("geolocation_data");
            v1 = geo_grp.getVar("latitude");
            v1.getVar(&latpix[0][0]);
            v1 = geo_grp.getVar("longitude");
            v1.getVar(&lonpix[0][0]);
            float fillval_geo, scale_factor;
            short senz_min, senz_max, sena_min, sena_max;
            NcVarAtt a1 = v1.getAtt("_FillValue");  // root group
            a1.getValues(&fillval_geo);
            v1 = geo_grp.getVar("sensor_azimuth_angle");
            v1.getVar(&senapix[0][0]);
            a1 = v1.getAtt("valid_min");  // root group
            a1.getValues(&sena_min);
            a1 = v1.getAtt("valid_max");  // root group
            a1.getValues(&sena_max);
            v1 = geo_grp.getVar("sensor_zenith_angle");
            v1.getVar(&senzpix[0][0]);
            a1 = v1.getAtt("valid_min");  // root group
            a1.getValues(&senz_min);
            a1 = v1.getAtt("valid_max");  // root group
            a1.getValues(&senz_max);
            a1 = v1.getAtt("scale_factor");  // root group
            a1.getValues(&scale_factor);

            senz_min *= scale_factor;
            senz_max *= scale_factor;
            sena_min *= scale_factor;
            sena_max *= scale_factor;

            if(l1cinput->verbose)cout << "computing parallax for each pixel.............." << endl;

            for (unsigned int i = 0; i < num_scans; i++) {
                for (unsigned int j = 0; j < num_pixels; j++) {
                    // Displacement vector -- wang et al. 2011 ----
                    temp = senapix[i][j] * scale_factor;
                    if (senzpix[i][j] * scale_factor >= senz_min &&
                        senzpix[i][j] * scale_factor <= senz_max && temp >= sena_min && temp <= sena_max) {
                        dv = Hsat * cth[i][j] * tan(senzpix[i][j] * scale_factor * M_PI / 180) /
                             (Hsat - cth[i][j]);
                        lat_new[i][j] = latpix[i][j] * M_PI / 180 - dv * cos(temp * M_PI / 180 + M_PI) / Re;
                        lon_new[i][j] =
                            lonpix[i][j] * M_PI / 180 -
                            dv * sin((temp * M_PI / 180 + M_PI)) / (Re * cos(lat_new[i][j] * M_PI / 180));
                        if (lon_new[i][j] * 180 / M_PI > 180) {
                            lon_new[i][j] -= 2 * M_PI;
                        }
                        if (lon_new[i][j] * 180 / M_PI < -180) {
                            lon_new[i][j] += 2 * M_PI;
                        }
                        lat_new[i][j] *= 180 / M_PI;
                        lon_new[i][j] *= 180 / M_PI;
                    } else {
                        lat_new[i][j] = fillval_geo;
                        lon_new[i][j] = fillval_geo;
                    }
                }
            }

            // write the new L1B with parallax correction
            string gran_str = to_string(gran);
            string cth_str = to_string(l1cinput->cloud_height);
            string l12_gran_new =
                l12_str.substr(0, 24) + "_CTH" + cth_str + ".L1B" + ".gran" + gran_str + ".nc";

            if(l1cinput->verbose)cout << "corrected L1B grid" << l12_gran_new << "     granule #" << gran_str << endl;

            const char* l12_gran = l12_gran_new.c_str();

            try {
                nc_l12_new = new NcFile(l12_gran, NcFile::replace);
            } catch (NcException& e) {
                e.what();
                cerr << "l1cgen mode= 9:: Failure write FULL L1B parallax corrected: " + l12_gran_new << endl;
                exit(1);
            }
            // global attributes--
            nc_l12_new->putAtt("title", "PACE OCIS Level-1B Data");
            nc_l12_new->putAtt("instrument", "OCIS");
            nc_l12_new->putAtt("processing_version", l1cinput->pversion);
            nc_l12_new->putAtt("Conventions", "CF-1.8 ACDD-1.3");
            nc_l12_new->putAtt("institution",
                               "NASA Goddard Space Flight Center, Ocean Biology Processing Group");
            nc_l12_new->putAtt(
                "license",
                "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/");
            nc_l12_new->putAtt("naming_authority", "gov.nasa.gsfc.sci.oceancolor");
            nc_l12_new->putAtt("keywords_vocabulary",
                               "NASA Global Change Master Directory (GCMD) Science Keywords");
            nc_l12_new->putAtt("stdname_vocabulary", "NetCDF Climate and Forecast (CF) Metadata Convention");
            nc_l12_new->putAtt("creator_name", "NASA/GSFC");
            nc_l12_new->putAtt("creator_email", "data@oceancolor.gsfc.nasa.gov");
            nc_l12_new->putAtt("creator_url", "https://oceancolor.gsfc.nasa.gov");
            nc_l12_new->putAtt("project", "PACE Project");
            nc_l12_new->putAtt("publisher_name", "NASA/GSFC");
            nc_l12_new->putAtt("publisher_email", "data@oceancolor.gsfc.nasa.gov");
            nc_l12_new->putAtt("publisher_url", "https://oceancolor.gsfc.nasa.gov");
            nc_l12_new->putAtt("processing_level", "L1B parallax");
            nc_l12_new->putAtt("cdm_data_type", "swath");
            nc_l12_new->putAtt("normalizedLt", "0LL");
            nc_l12_new->putAtt("history", l1cinput->history);
            nc_l12_new->putAtt("product_name", l12_gran_new);
            NcGroupAtt i1 = nc_l12->getAtt("startDirection");  // NcGroupAtt is a global attr!!
            i1.getValues(name);
            nc_l12_new->putAtt("startDirection", name);
            i1 = nc_l12->getAtt("endDirection");
            i1.getValues(name);
            nc_l12_new->putAtt("endDirection", name);
            nc_l12_new->putAtt("sample_offset", "0LL");
            i1 = nc_l12->getAtt("time_coverage_start");
            string timeStartStr;
            i1.getValues(timeStartStr);
            nc_l12_new->putAtt("time_coverage_start", timeStartStr);
            i1 = nc_l12->getAtt("time_coverage_end");
            i1.getValues(name);
            nc_l12_new->putAtt("time_coverage_end", name);
            i1 = nc_l12->getAtt("date_created");
            i1.getValues(name);
            nc_l12_new->putAtt("date_created", name);
            nc_l12_new->putAtt("earth_sun_distance_correction",ncFloat,FillValue);

            NcDim ydim = nc_l12_new->addDim("number_of_scans", num_scans);
            NcDim xdim = nc_l12_new->addDim("ccd_pixels", num_pixels);

            std::vector<NcDim> dimvec2_geo;
            dimvec2_geo.push_back(ydim);
            dimvec2_geo.push_back(xdim);

            nc_l12_new->addDim("SWIR_pixels", num_pixels);
            num_blue_bands = 120;
            num_red_bands = 120;
            num_SWIR_bands = 9;
            size_t vector_elements = 3;
            nc_l12_new->addDim("blue_bands", num_blue_bands);
            nc_l12_new->addDim("red_bands", num_red_bands);
            nc_l12_new->addDim("SWIR_bands", num_SWIR_bands);
            nc_l12_new->addDim("vector_elements", vector_elements);

            ba_grp = nc_l12_new->addGroup("bin_attributes");
            geo_grp = nc_l12_new->addGroup("geolocation_data");
            v1 = ba_grp.addVar("nadir_view_time", ncDouble, ydim);

            string longName = "Time bin was viewed at nadir view";
            v1.putAtt("long_name", longName);
            string units = "seconds since " + timeStartStr.substr(0, 10); // just add the date, no time
            v1.putAtt("units", units);
            v1.putAtt("_FillValue", ncDouble, FillValue2);
            double valid_min_d = 0;
            v1.putAtt("valid_min", ncDouble, valid_min_d);
            double valid_max_d = 172800;
            v1.putAtt("valid_max", ncDouble, valid_max_d);
            v1.putVar(&ptime[0]);

            v1 = geo_grp.addVar("latitude", ncFloat, dimvec2_geo);
            longName = "Latitudes of bin locations";
            v1.putAtt("long_name", longName);
            units = "degrees_north";
            v1.putAtt("units", units);
            v1.putAtt("_FillValue", ncFloat, FillValue);
            float valid_min = -90;
            v1.putAtt("valid_min", ncFloat, valid_min);
            float valid_max = 90;
            v1.putAtt("valid_max", ncFloat, valid_max);
            v1.putVar(&lat_new[0][0]);

            v1 = geo_grp.addVar("longitude", ncFloat, dimvec2_geo);
            longName = "Longitudes of bin locations";
            v1.putAtt("long_name", longName);
            units = "degrees_east";
            v1.putAtt("units", units);
            v1.putAtt("_FillValue", ncFloat, FillValue);
            valid_min = -180;
            v1.putAtt("valid_min", ncFloat, valid_min);
            valid_max = 180;
            v1.putAtt("valid_max", ncFloat, valid_max);
            v1.putVar(&lon_new[0][0]);

            nc_l12->close();
            nc_l12_new->close();

            free (lat_new);
            free (lon_new);
            free (ptime);
            free (latpix);
            free (lonpix);
            free (senapix);
            free (senzpix);
            free (cth);
        }
    } else {
        cout << "l1cinput->cloud_correct WAS NOT setup!" << endl;
        exit(1);
    }

    return 0;
}



int32_t L1C::swtime_swt2_segment(int swt, L1C_input* l1cinput, l1c_filehandle* l1cfile, int32_t norbs,
                                 double* tswt, double tcross, double mgv, double* tmgv,double *orb_time_tot,size_t norbs_tot) {
    int16_t bina = 0, binb = 0, ngridlines, gd = 0;
    double tg = 0., mot = 0, t_start = -1, t_end = -1;
    int16_t gn=-1,gn2=-1;
    t_start = tswt[0];
    t_end = tswt[norbs - 1];

    mot = ((l1cinput->grid_resolution) * 1000) / mgv;  // in seconds
    l1cfile->mot=mot;
    if(tcross>0)
    {
      if(orb_time_tot[0]>tcross)
      { 
      gn=(tcross+24*3600-orb_time_tot[0])/l1cfile->mot;
      gn2=(orb_time_tot[norbs_tot-1]-tcross)/l1cfile->mot;
      }
      else if(orb_time_tot[norbs_tot-1]<tcross)
      {
      gn=(tcross-orb_time_tot[0])/l1cfile->mot;
      gn2=(orb_time_tot[norbs_tot-1]+24*3600-tcross)/l1cfile->mot;
      }
      else
      {
      gn=(tcross-orb_time_tot[0])/l1cfile->mot;
      gn2=(orb_time_tot[norbs_tot-1]-tcross)/l1cfile->mot;
      }

     if(l1cinput->verbose)
     {
     cout<<"tcross "<<tcross<<"orb_time_tot[ini] "<<orb_time_tot[0]<<"gn "<<gn<<"mot "<<l1cfile->mot<<"gn2 "<<gn2<<endl;
     cout<<"#gridlines "<<gn+gn2<<"orb_time_tot[end] "<<orb_time_tot[norbs_tot-1]<<endl;
     }


    l1cfile->num_gridlines=gn+gn2;
   }

    
    ngridlines = l1cfile->num_gridlines;

    int flag_time = -1;

    if (tcross > 0.0) {
        flag_time = 0;
        if(l1cinput->verbose)cout << "computing time series assuming mean swath velocity..for swath#." << swt << endl;

        // backward section of swath--before equat
        tg = tcross - mot / 2;
        gd = ngridlines / 2 - 1;

        tmgv[gd] = tg;
        binb = 1;

        while (tg >= t_start) {

            binb++;
            tg -= mot;
            gd--;
            tmgv[gd] = tg;

        }
        // forward section of swath--after equator
        tg = tcross + mot / 2;
        gd = ngridlines / 2;

        tmgv[gd] = tg;
        bina = 1;

        while (tg <= t_end) {

            bina++;
            tg += mot;
            gd++;
            tmgv[gd] = tg;

        }

        if(l1cinput->verbose)cout << "bina.." << bina << "binb.." << binb << endl;

        l1cfile->num_gridlines = binb + bina;

        if(l1cinput->verbose)cout << "number of L1C gridlines along-track..." << l1cfile->num_gridlines << "for swath #.." << swt
             << "t_start.." << t_start << "t_end..." << t_end << endl;
    }  // end equat crossing
    else {
        if(l1cinput->verbose)cout << "time series not possible for swath #.." << swt << "tcross<0...." << endl;
        flag_time = 1;
    }

    return flag_time;
}

int32_t L1C::swtime_swt2(int swt, L1C_input* l1cinput, l1c_filehandle* l1cfile, int32_t norbs, double* tswt,
                         double tcross, double mgv, double* tmgv, double *orb_time_tot, size_t norbs_tot) {
    int16_t bina = 0, binb = 0, gd = 0;
    double tg = 0., mot = 0;
    int16_t gn=-1,gn2=-1;
    mot = ((l1cinput->grid_resolution) * 1000) / mgv;  // in seconds
    l1cfile->mot=mot;
    if(tcross>0)
    {
    if(orb_time_tot[0]>tcross)
      {
      gn=(tcross+24*3600-orb_time_tot[0])/l1cfile->mot;
      gn2=(orb_time_tot[norbs_tot-1]-tcross)/l1cfile->mot;
      }
      else if(orb_time_tot[norbs_tot-1]<tcross)
      {
      gn=(tcross-orb_time_tot[0])/l1cfile->mot;
      gn2=(orb_time_tot[norbs_tot-1]+24*3600-tcross)/l1cfile->mot;
      }
      else
      {
      gn=(tcross-orb_time_tot[0])/l1cfile->mot;
      gn2=(orb_time_tot[norbs_tot-1]-tcross)/l1cfile->mot;
      }

       if(l1cinput->verbose)
       {cout<<"tcross "<<tcross<<"orb_time_tot[ini] "<<orb_time_tot[0]<<"gn "<<gn<<"mot "<<l1cfile->mot<<"gn2 "<<gn2<<endl;
                cout<<"#gridlines "<<gn+gn2<<"orb_time_tot[end] "<<orb_time_tot[norbs_tot-1]<<endl;
       }

    l1cfile->num_gridlines=gn+gn2;  
    }




    int flag_time = -1;

    if (tcross > 0.0) {
        flag_time = 0;
        if(l1cinput->verbose)cout << "computing time series assuming mean swath velocity..for swath#." << swt << endl;

        tg = tcross - mot / 2;

        gd=gn-1;
        tmgv[gd] = tg;
        binb = 1;

        if(tg<0 && l1cinput->verbose){
            cout<<"NEGATIVE FIRST TIME BEFORE CROSSING"<<tg<<"in swtime_swt2 for mean vel tmgv calc "<<mgv<<"mot/2 "<<mot/2<<"tcross "<<tcross<<endl;
        }

        while (binb < gn) {
            binb++;
            tg -= mot;
            gd--;
            tmgv[gd] = tg;

        }
        tg = tcross + mot / 2;

        gd=gn;
        tmgv[gd] = tg;
        bina = 1;


        while (bina < gn2+1) {  
            bina++;
            tg += mot;
            gd++;
            tmgv[gd] = tg;

        }

        if(l1cinput->verbose)cout << "bina.." << bina << "binb.." << binb << endl;

        l1cfile->num_gridlines = binb + bina;

        if(l1cinput->verbose)cout << "number of L1C gridlines along-track..." << l1cfile->num_gridlines << "for swath #.." << swt
             << endl;
    }  // end equat crossing
    else {
        if(l1cinput->verbose)cout << "time series not possible for swath #.." << swt << "tcross<0...." << endl;
        flag_time = 1;
    }

    return flag_time;
}


int32_t L1C::ect_swt(l1c_filehandle* l1cfile, int ix1, int ix2, double* tswt_tot, double* latswt_tot,
                     double* lonswt_tot, double* ovel_tot, double* gvel_tot, double* tswt, double* latswt,
                     double* lonswt, float* tcross, float* loncross, double* ovel, double* gvel) {
    // determine equatorial crossing time--------------
    double t1 = -1, t2 = -1;
    int asc_mode = -1;
    size_t tindex = -1;
    int scan_brack[2] = {-1, -1};
    float lat1, lat2, lon1, lon2;  // lat,lon defined globally
    int flag_ect = -1;
    double sum1 = 0, sum3 = 0;
    size_t c = 0, norbs;

    norbs = ix2 - ix1 + 1;
    int kk = ix1;

    while (kk <= ix2) {
        latswt[c] = latswt_tot[kk];
        lonswt[c] = lonswt_tot[kk];
        tswt[c] = tswt_tot[kk];

        sum1 += ovel_tot[kk];
        sum3 += gvel_tot[kk];
        kk++;
        c++;
    }

    // determine orbit direction--asc/desc
    if (latswt[0] > latswt[1])
        asc_mode = 0;  // descending
    else
        asc_mode = 1;  // ascending


    for (size_t i = 0; i < norbs - 1; i++) {  // improve with binary search
        if (latswt[i] < 0. && latswt[i + 1] > 0.) {//ascending orbit
            asc_mode=1;
            scan_brack[0] = i + 1;
            scan_brack[1] = i + 2;
            lat1 = latswt[i];
            lat2 = latswt[i + 1];
            lon1 = lonswt[i];
            lon2 = lonswt[i + 1];
            if(i<norbs - 2 && i>0)
            {
            if(lat1==latswt[i-1] || lat2==latswt[i + 2])
            scan_brack[0] = -1;
            scan_brack[1] = -1;
            }
            
            tindex = i;
            i = norbs;           
            break;
        }

        else if (latswt[i] > 0. && latswt[i + 1] < 0) {  // descending orbit
            asc_mode=0;
            scan_brack[0] = i + 2;                                        // negative lat first convention
            scan_brack[1] = i + 1;
            lat1 = latswt[i + 1];
            lat2 = latswt[i];
            lon1 = lonswt[i + 1];
            lon2 = lonswt[i];
            if(i<norbs - 2 && i>0)
            {          
            if(lat2==latswt[i-1] || lat1==latswt[i + 2])
            scan_brack[0] = -1;
            scan_brack[1] = -1;
            }

            tindex = i;
            i = norbs;
            break;
        }
    }  // end for

    if(l1cfile->verbose)        cout<<"computing equatorial crossing time and longitude at crossing for a specific swath..orbit direction 0: descending, 1: ascending..."<<asc_mode<<endl;
    // interpolate time--linear
    if (scan_brack[0] > -1) {  // if file with equat crossing
        t1 = tswt[tindex];
        t2 = tswt[tindex + 1];
        if(l1cfile->verbose)cout << "lat1.." << lat1 << "lat2..." << lat2 << "t1.." << t1 << "t2..." << t2 << endl;
        *tcross = t1 - lat1 * (t2 - t1) / (lat2 - lat1);  // equat crossing time

        float dtcross = (*tcross - t1) / (t2 - t1);
        *loncross = lon1 + (lon2 - lon1) * dtcross;  // latcross is zero

        l1cfile->eqt = *tcross;
        l1cfile->orbit_node_lon = *loncross;
        l1cfile->orb_dir = asc_mode;
        flag_ect = 0;
    } else {
        l1cfile->eqt = -999.0;
        l1cfile->orbit_node_lon = -999.0;
        l1cfile->orb_dir = asc_mode;
        flag_ect = 1;
    }

    // computing mean swath velocity (m/s) --------------------------

    *ovel = (sum1 / norbs);         // mean velocity in m/s
    *gvel = (sum3 / norbs) * 1000;  // in meters/s

    return flag_ect;
}

int sunz_swt(int ix1_swt, int ix2_swt, int16_t* hour_swt, int16_t* day_swt, int16_t* year_swt, double* olat,
             double* olon) {
    int day_mode = -1;  // daytime is 1, nighttime is 0
    float sunz1, suna1;
    int year1 = year_swt[ix1_swt];
    int day1 = day_swt[ix1_swt];
    float hour1 = hour_swt[ix1_swt];
    float lat1 = olat[ix1_swt];
    float lon1 = olon[ix1_swt];

    sunangs_(&year1, &day1, &hour1, &lon1, &lat1, &sunz1, &suna1);
    if (sunz1 <= 93)
        day_mode = 1;
    else
        day_mode = 0;  // nighttime

 //   if(l1cfile->verbose)cout << "lat1.." << lat1 << "lon1.." << lon1 << "year1" << year1 << "day1" << day1 << "hour1" << hour1<< "sunz1..first index" << sunz1 << "ix1" << ix1_swt << "ix2" << ix2_swt << endl;

    year1 = year_swt[ix2_swt];
    day1 = day_swt[ix2_swt];
    hour1 = hour_swt[ix2_swt];
    lat1 = olat[ix2_swt];
    lon1 = olon[ix2_swt];

    sunangs_(&year1, &day1, &hour1, &lon1, &lat1, &sunz1, &suna1);
    if (sunz1 <= 93 || day_mode == 1)
        day_mode = 1;
    else
        day_mode = 0;

  // if(l1cfile->verbose)cout << "lat1.." << lat1 << "lon1.." << lon1 << "year1" << year1 << "day1" << day1 << "hour1" << hour1<< "sunz1..second index" << sunz1 << "ix1" << ix1_swt << "ix2" << ix2_swt << endl;

    return day_mode;
}

int32_t L1C::create_time_swt(int num_gridlines, double tfile_ini_sec, double* tmgvf, double tswt_ini_sec,
                             double tswt_end_sec, string* tswt_ini, string* tswt_ini_file, string* tswt_mid,string* tswt_end) {
    int16_t oneyear, onemon, oneday, onehour, onemin;
    double onesec,tswt_mid_sec;
    string y_swt, mo_swt, d_swt, h_swt, mi_swt, s_swt, s_swt2;
    int logoff = -1;

    unix2ymdhms(tswt_ini_sec, &oneyear, &onemon, &oneday, &onehour, &onemin, &onesec);
    y_swt = std::to_string(oneyear);
    mo_swt = std::to_string(onemon);
    d_swt = std::to_string(oneday);
    h_swt = std::to_string(onehour);
    mi_swt = std::to_string(onemin);
    s_swt2 = std::to_string(round(onesec));

    length = (int)floor(log10(onemon)) + 1;
    if (length == 1)
        mo_swt = "0" + mo_swt;
    length = (int)floor(log10(oneday)) + 1;
    if (length == 1)
        d_swt = "0" + d_swt;
    if (onehour == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(onehour + logoff)) + 1;
    if (length == 1)
        h_swt = "0" + h_swt;
    if (onemin == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(onemin + logoff)) + 1;
    if (length == 1)
        mi_swt = "0" + mi_swt;
    if (onesec == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(round(onesec) + logoff)) + 1;
    if (length == 1)
        s_swt2 = "0" + s_swt2;
    if (s_swt2.substr(1, 1) == ".")
        s_swt2 = "00";

    *tswt_ini = y_swt + "-" + mo_swt + "-" + d_swt + "T" + h_swt + ":" + mi_swt + ":" + s_swt2.substr(0, 2);
    *tswt_ini_file = y_swt + mo_swt + d_swt + "T" + h_swt + mi_swt + s_swt2.substr(0, 2);

    tswt_end_sec = tfile_ini_sec + tmgvf[num_gridlines - 2];
    unix2ymdhms(tswt_end_sec, &oneyear, &onemon, &oneday, &onehour, &onemin, &onesec);

    y_swt = std::to_string(oneyear);
    mo_swt = std::to_string(onemon);
    d_swt = std::to_string(oneday);
    h_swt = std::to_string(onehour);
    mi_swt = std::to_string(onemin);
    s_swt2 = std::to_string(round(onesec));

    length = (int)floor(log10(onemon)) + 1;
    if (length == 1)
        mo_swt = "0" + mo_swt;
    length = (int)floor(log10(oneday)) + 1;
    if (length == 1)
        d_swt = "0" + d_swt;
    if (onehour == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(onehour + logoff)) + 1;
    if (length == 1)
        h_swt = "0" + h_swt;
    if (onemin == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(onemin + logoff)) + 1;
    if (length == 1)
        mi_swt = "0" + mi_swt;
    if (onesec == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(round(onesec) + logoff)) + 1;
    if (length == 1)
        s_swt2 = "0" + s_swt2;
    if (s_swt2.substr(1, 1) == ".")
        s_swt2 = "00";

    *tswt_end = y_swt + "-" + mo_swt + "-" + d_swt + "T" + h_swt + ":" + mi_swt + ":" + s_swt2.substr(0, 2);

    tswt_mid_sec=(tswt_ini_sec+tswt_end_sec)/2;
    unix2ymdhms(tswt_mid_sec, &oneyear, &onemon, &oneday, &onehour, &onemin, &onesec);

    y_swt = std::to_string(oneyear);
    mo_swt = std::to_string(onemon);
    d_swt = std::to_string(oneday);
    h_swt = std::to_string(onehour);
    mi_swt = std::to_string(onemin);
    s_swt2 = std::to_string(round(onesec));

    length = (int)floor(log10(onemon)) + 1;
    if (length == 1)
        mo_swt = "0" + mo_swt;
    length = (int)floor(log10(oneday)) + 1;
    if (length == 1)
        d_swt = "0" + d_swt;
    if (onehour == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(onehour + logoff)) + 1;
    if (length == 1)
        h_swt = "0" + h_swt;
    if (onemin == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(onemin + logoff)) + 1;
    if (length == 1)
        mi_swt = "0" + mi_swt;
    if (onesec == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(round(onesec) + logoff)) + 1;
    if (length == 1)
        s_swt2 = "0" + s_swt2;
    if (s_swt2.substr(1, 1) == ".")
        s_swt2 = "00";

    *tswt_mid = y_swt + "-" + mo_swt + "-" + d_swt + "T" + h_swt + ":" + mi_swt + ":" + s_swt2.substr(0, 2);

    return 0;
}

// open telemetry file parameters for creating L1C grid----
int32_t L1C::open_l1atol1c3(L1C_input* l1cinput, l1c_filehandle* l1cfile) {
    int32_t n_files;
    const char* ptstr;
    char* ifile_char;
    string str, ifile_str;
    int status = -1, status1 = -1, status2 = -1;
    unsigned char** hkpackets = nullptr;
    uint8_t* apids = nullptr;
    size_t number_hkpack, number_scpack, number_orecords, number_hkpack_tot=0, number_scpack_tot=0,
        number_orecords_tot, nr=0;
    uint8_t ubnum1, ubnum2;
    double *sec = nullptr, *tai58_sec = nullptr;
    int16_t *year = nullptr, *mon = nullptr, *day = nullptr, *hour = nullptr, *min = nullptr;
    int16_t oneyear, onemon, oneday, onehour, onemin;
    double onesec;
    double *orb_time = nullptr, *orb_lat = nullptr, *orb_lon = nullptr;
    float **lat_gd = nullptr, **lon_gd = nullptr, **alt = nullptr;
    double omeg = 7.29211585494e-5;
    short n_ephem;
    string temp_str, tai_str;
    double vxyz = 0, mov1 = 0., mgv1 = 0., mov2 = 0, mgv2, *orb_vel = nullptr, *grn_vel = nullptr;
    int* orb_dir = nullptr;
    int ix1 = -1, ix2 = -1, ix3 = -1, ix4 = -1, ix5 = -1, ix6 = -1;
    double *tmgv1 = nullptr, *tmgv2 = nullptr, *tmgvf = nullptr, *tmgvf2 = nullptr;
    int32_t num_gridlines = -1, norbs = -1;
    double *tswt = nullptr, *latswt = nullptr, *lonswt = nullptr;
    double *tswt2 = nullptr, *latswt2 = nullptr, *lonswt2 = nullptr;
    const char* outlist;
    string ofile_str, senstr;
    string y_swt, mo_swt, d_swt, h_swt, mi_swt, s_swt, tswt_ini, tswt_end,tswt_mid;
    int32_t gd_per_gran = -1, numgran = -1;
    double deltasec = -1;
    string s_swt2;
    string tswt_ini_file;
    double rl2, pos_norm, clatg2, fe = 1 / 298.257;
    int day_mode = -1;  // 0 is nighttime 1 is dayttime
    double tai58unix;
    double tswt_ini_sec, tswt_end_sec;
    vector<size_t> ix;
    double tfile_ini, tfile_end;
    int16_t syear, smon, sday, shour, smin, syear2, smon2, sday2, shour2, smin2;
    double secs, second, secs2, second2;
    double tfile_ini_sec, tfile_end_sec;
    size_t ix_ini=-1,ix_end=-1;

    // TOTAL ARRAYS---ASSUMING 6000 rcords
    number_orecords_tot = 6000;
    size_t cc = 0, c = 0;

    int16_t* year_tot = (int16_t*)calloc(number_orecords_tot, sizeof(int16_t));
    int16_t* day_tot = (int16_t*)calloc(number_orecords_tot, sizeof(int16_t));
    int16_t* hour_tot = (int16_t*)calloc(number_orecords_tot, sizeof(int16_t));
    double* orb_time_tot = (double*)calloc(number_orecords_tot, sizeof(double));
    double* orb_lat_tot = (double*)calloc(number_orecords_tot, sizeof(double));
    double* orb_lon_tot = (double*)calloc(number_orecords_tot, sizeof(double));
    int* orb_dir_tot = (int*)calloc(number_orecords_tot, sizeof(int));
    orb_array2* posr_tot = new orb_array2[number_orecords_tot]();
    orb_array2* velr_tot = new orb_array2[number_orecords_tot]();
    double* orb_vel_tot = (double*)calloc(number_orecords_tot, sizeof(double));
    double* grn_vel_tot = (double*)calloc(number_orecords_tot, sizeof(double));

    //*************************************************************
    ifile_str = l1cinput->files[0];
    ifile_char = (char*)ifile_str.c_str();
    ofile_str = ifile_str.substr(0, 24);
    outlist = "list_l1c_granules.txt";  // DEFAULT OFILE

    if (l1cinput->outlist[0] == '\0') {  // only when no txt output list file is provided
        strcpy(l1cinput->outlist, outlist);
        if(l1cinput->verbose)cout << "L1C granules written to DEFAULT file........" << outlist << endl;
    } else {
        if(l1cinput->verbose)cout << "L1C granules written to "
                "file......................................................................................."
             << l1cinput->outlist << endl;
    }

    l1cfile->gransize = l1cinput->gransize;

    file_format format = getFormat(ifile_char);
    if(l1cinput->verbose)cout << "format.type.." << format.type << endl;

    l1cfile->sensor = l1cinput->sensor;  // SPEX 1, OCI 2 and HARP 3

    if (l1cinput->sensor == 34) {
        senstr = "SPEXONE";
        l1cfile->nbinx = 29;
        l1cfile->n_views = 10;
    } else if (l1cinput->sensor == 30) {
        senstr = "OCI";
        l1cfile->nbinx = 519;
        l1cfile->n_views = 2;
    } else if (l1cinput->sensor == 35) {
        senstr = "HARP2";
        l1cfile->nbinx = 457;
        l1cfile->n_views = 90;
    } else {
        if(l1cinput->verbose)cout << "sensor by default is OCI option 2....." << endl;
        senstr = "OCI";
        l1cfile->nbinx = 514;
        l1cfile->n_views = 2;
    }

    if(l1cinput->verbose)cout << "PACE sensor to be griddded....." << senstr << endl;

    n_files = l1cfile->ifiles.size();
    if(l1cinput->verbose)cout << "number of files in the list...." << n_files << endl;

    int32_t nfiles = l1cfile->ifiles.size();
  
    int FirsTerrain=0;

    for (int fi = 0; fi < nfiles; fi++) 
    {
        if(fi==nfiles-1)FirsTerrain=1;
        str = l1cfile->ifiles[fi];
        ptstr = str.c_str();
        if(l1cinput->verbose)
        {
        cout << "********************************************************************************************"
                "*************"
             << endl;
        cout << "Opening L1A file ..." << ptstr << "for telemetry......." << endl;
        cout << "********************************************************************************************"
                "*************"
             << endl;
        }

        NcFile* nc_l1a;
        try {
            nc_l1a = new NcFile(ptstr, NcFile::read);
        } catch (NcException& e) {
            e.what();
            cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: " + str << endl;
            exit(1);
        }

        string name;
        NcGroupAtt i1 = nc_l1a->getAtt("time_coverage_start");  // NcGroupAtt is a global attr!!
        i1.getValues(name);

        tfile_ini = isodate2unix(name.c_str());
       
        i1 = nc_l1a->getAtt("time_coverage_end");  // NcGroupAtt is a global attr!!
        i1.getValues(name);
        tfile_end = isodate2unix(name.c_str());

        // this is default start time, start_timeflag=1 or using start time of swath file HKT------
        unix2ymds(tfile_ini, &syear, &smon, &sday, &secs);
        unix2ymds(tfile_end, &syear2, &smon2, &sday2, &secs2);

        if(l1cinput->verbose)cout << "secs elapsed.." << secs << "initial granule #..." << round(secs / (l1cfile->gransize * 60))
             << endl;
        unix2ymdhms(tfile_ini, &syear, &smon, &sday, &shour, &smin, &second);
        if(l1cinput->verbose)cout << "HKT file start time................."
             << "year.." << syear << "month..." << smon << "day..." << sday << "hour.." << shour << "min.."
             << smin << "sec..." << second << endl;
        unix2ymdhms(tfile_end, &syear2, &smon2, &sday2, &shour2, &smin2, &second2);
        if(l1cinput->verbose)cout << "HKT file end time................."
             << "year.." << syear2 << "month..." << smon2 << "day..." << sday2 << "hour.." << shour2
             << "min.." << smin2 << "sec..." << second2 << endl;

        if (l1cinput->start_timeflag == 0) {  // forcing to hms =0:0:0
            secs = 0.;
           }
     
         tfile_ini_sec = ymds2unix(syear, smon, sday, secs);
         tfile_end_sec =
            tfile_ini_sec + 49 * 60 * 3;  // 49 minutes per half-orbit, 3x cause 1.5 orbit, forthy swaths (3 before)
        
        int16_t y_ini, mo_ini, d_ini, h_ini, mi_ini;
        double sec_ini;
        unix2ymdhms(tfile_ini_sec, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
        if(l1cinput->verbose)cout << "tfile_ini_sec.."
             << "YEAR.." << y_ini << "MONTH.." << mo_ini << "DAY.." << d_ini << "HOUR.." << h_ini << "MIN.."
             << mi_ini << "SEC.." << sec_ini << endl;
        unix2ymdhms(tfile_end_sec, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
        if(l1cinput->verbose)cout << "tfile_end_sec.."
             << "YEAR.." << y_ini << "MONTH.." << mo_ini << "DAY.." << d_ini << "HOUR.." << h_ini << "MIN.."
             << mi_ini << "SEC.." << sec_ini << endl;

        l1cfile->tfile_ini_sec = tfile_ini_sec;
        l1cfile->tfile_end_sec = tfile_end_sec;
        
        // open dimensions
        /// number of packets
        NcDim dimhp = nc_l1a->getDim("SC_hkt_pkts");
        NcDim dimtp = nc_l1a->getDim("max_SC_packet");
        NcDim dimor = nc_l1a->getDim("orb_records");

        number_hkpack = dimhp.getSize();
        number_scpack = dimtp.getSize();
        number_orecords = dimor.getSize();

        // total counters---
        number_hkpack_tot += number_hkpack;
        number_scpack_tot += number_scpack;
        nr += number_orecords;

        if(l1cinput->verbose)cout << "number_hkpack_tot.." << number_hkpack_tot << "number_scpack_tot.." << number_scpack_tot
             << "number of orbit records total.REAL." << nr << endl;

        // allocat mem
        hkpackets =
            allocate2d_uchar(number_hkpack, number_scpack);  // NUMBER of ephem elements x max SIZE of package
        apids = (uint8_t*)calloc(number_hkpack, sizeof(uint8_t));

        orb_time = (double*)calloc(number_orecords, sizeof(double));
        orb_lat = (double*)calloc(number_orecords, sizeof(double));
        orb_lon = (double*)calloc(number_orecords, sizeof(double));
        orb_dir = (int*)calloc(number_orecords, sizeof(int));

        // open groups
        NcGroup telGrp = nc_l1a->getGroup("housekeeping_data");
        NcGroup navGrp = nc_l1a->getGroup("navigation_data");

        // open vars ids
        NcVar v1 = telGrp.getVar("SC_HKT_packets");
        v1.getVar(&hkpackets[0][0]);

        v1 = navGrp.getVar("orb_time");
        v1.getVar(&orb_time[0]);
        v1 = navGrp.getVar("orb_lat");
        v1.getVar(&orb_lat[0]);
        v1 = navGrp.getVar("orb_lon");
        v1.getVar(&orb_lon[0]);

        for (size_t hk = 0; hk < number_hkpack; hk++) {
            ubnum1 = (uint8_t)hkpackets[hk][0];  //-48;//48 is 0 ascii code
            ubnum2 = (uint8_t)hkpackets[hk][1];

            apids[hk] = (ubnum1 % 8) * 256 + ubnum2;
            if (apids[hk] == 128) {
                ix.push_back(hk);  // packet index where ephemr
            }
        }

        if(l1cinput->verbose)cout << "#number of ephem elements...." << ix.size() << "for HKT file..." << ptstr << "#orecords..."
             << number_orecords << endl;
        n_ephem = ix.size();

        // Reverse the byte order  from small endian to large endian, small endian last byte is stored first
        // convert to double and later swap endian

        // allocate mem
        orb_vel = (double*)calloc(n_ephem, sizeof(double));
        grn_vel = (double*)calloc(n_ephem, sizeof(double));

        sec = (double*)calloc(n_ephem, sizeof(double));
        year = (int16_t*)calloc(n_ephem, sizeof(int16_t));
        mon = (int16_t*)calloc(n_ephem, sizeof(int16_t));
        day = (int16_t*)calloc(n_ephem, sizeof(int16_t));
        hour = (int16_t*)calloc(n_ephem, sizeof(int16_t));
        min = (int16_t*)calloc(n_ephem, sizeof(int16_t));

        tai58_sec = (double*)calloc(n_ephem, sizeof(double));

        // process all packets----
        double tai58;
        c = 0;
        // TIME

        while (c < ix.size()) {
            double* tai_ptr = (double*)(hkpackets[ix[c]] + 16);  // moving pointer to position 16
            swapc_bytes((char*)tai_ptr, 8,
                        1);  // swap 8 bytes once or 1 16-23 pos, lets say we have 16-32 indexes and we want
                             // to swap 8 bytes at the time, some we need ntime=2
            tai58 = *tai_ptr;
            tai58unix = tai58_to_unix(tai58);

            unix2ymdhms(tai58unix, &oneyear, &onemon, &oneday, &onehour, &onemin, &onesec);

            year[c] = oneyear;
            mon[c] = onemon;
            day[c] = oneday;
            hour[c] = onehour;
            min[c] = onemin;
            sec[c] = onesec;

            tai58_sec[c] = tai58unix;

            c++;
        }

        // POS/VEL alloc mem-----
        orb_array2* posr = new orb_array2[n_ephem]();
        orb_array2* velr = new orb_array2[n_ephem]();

        double dp1, dp2, dp3;

        // computed orbital velocity from ephemeris
        c = 0;
        while (c < ix.size()) {
            double* posi = (double*)(hkpackets[ix[c]] + 120);
            double* veli = (double*)(hkpackets[ix[c]] + 144);
            double* ecmat = (double*)(hkpackets[ix[c]] + 176);

            swapc_bytes((char*)posi, 8, 3);
            swapc_bytes((char*)veli, 8, 3);
            swapc_bytes((char*)ecmat, 8, 9);

            // assuming ecmat index 1-3 first row, 4-6 second row and 7-9 third row
            // dot product
            // position
            dp1 = posi[0] * ecmat[0] + posi[1] * ecmat[1] + posi[2] * ecmat[2];
            dp2 = posi[0] * ecmat[3] + posi[1] * ecmat[4] + posi[2] * ecmat[5];
            dp3 = posi[0] * ecmat[6] + posi[1] * ecmat[7] + posi[2] * ecmat[8];

            posr[c][0] = dp1;
            posr[c][1] = dp2;
            posr[c][2] = dp3;

            // velocity
            dp1 = veli[0] * ecmat[0] + veli[1] * ecmat[1] + veli[2] * ecmat[2];
            dp2 = veli[0] * ecmat[3] + veli[1] * ecmat[4] + veli[2] * ecmat[5];
            dp3 = veli[0] * ecmat[6] + veli[1] * ecmat[7] + veli[2] * ecmat[8];

            velr[c][0] = dp1;
            velr[c][1] = dp2;
            velr[c][2] = dp3;

            velr[c][0] = velr[c][0] + posr[c][1] * omeg;
            velr[c][1] = velr[c][1] - posr[c][0] * omeg;

            // computing orbital velocity
            vxyz = sqrt(velr[c][0] * velr[c][0] + velr[c][1] * velr[c][1] +
                        velr[c][2] * velr[c][2]);  // units m/s
            orb_vel[c] = vxyz;

            pos_norm = sqrt(posr[c][0] * posr[c][0] + posr[c][1] * posr[c][1] + posr[c][2] * posr[c][2]);
            clatg2 = sqrt(posr[c][0] * posr[c][0] + posr[c][1] * posr[c][1]) / pos_norm;
            rl2 = Re * (1 - fe) / (sqrt(1 - (2 - fe) * fe * clatg2 * clatg2));
            grn_vel[c] = vxyz * rl2 / pos_norm;

            c++;
        }

        for (short n = 0; n < n_ephem; n++) {
            orb_time_tot[cc] = orb_time[n];
            orb_lon_tot[cc] = orb_lon[n];
            orb_lat_tot[cc] = orb_lat[n];
            orb_vel_tot[cc] = orb_vel[n];
            grn_vel_tot[cc] = grn_vel[n];
            velr_tot[cc][0] = velr[n][0];
            velr_tot[cc][1] = velr[n][1];
            velr_tot[cc][2] = velr[n][2];
            posr_tot[cc][0] = posr[n][0];
            posr_tot[cc][1] = posr[n][1];
            posr_tot[cc][2] = posr[n][2];
            year_tot[cc] = year[n];
            day_tot[cc] = day[n];
            hour_tot[cc] = hour[n];        
            cc++;
        }

        if (hkpackets != nullptr)
            free (hkpackets);
        hkpackets = nullptr;

        if (apids != nullptr)
            free (apids);
        apids = nullptr;

        ix.clear();

        if (tai58_sec != nullptr)
            free (tai58_sec);
        tai58_sec = nullptr;

        if (sec != nullptr)
            free (sec);
        sec = nullptr;

        if (year != nullptr)
            free (year);
        year = nullptr;

        if (mon != nullptr)
            free (mon);
        mon = nullptr;

        if (day != nullptr)
            free (day);
        day = nullptr;

        if (hour != nullptr)
            free (hour);
        hour = nullptr;

        if (min != nullptr)
            free (min);
        min = nullptr;

        if (orb_lat != nullptr)
            free (orb_lat);
        orb_lat = nullptr;

        if (orb_lon != nullptr)
            free (orb_lon);
        orb_lon = nullptr;

        if (orb_time != nullptr)
            free (orb_time);
        orb_time = nullptr;

        if (orb_vel != nullptr)
            free (orb_vel);
        orb_vel = nullptr;

        if (grn_vel != nullptr)
            free (grn_vel);
        grn_vel = nullptr;

        if (orb_dir != nullptr)
            free (orb_dir);
        orb_dir = nullptr;

        if (posr != nullptr)
            delete[] (posr);
        posr = nullptr;

        if (velr != nullptr)
            delete[] (velr);
        velr = nullptr;
        if (senz != nullptr)
            free (senz);
        senz = nullptr;

        if (sena != nullptr)
            free (sena);
        sena = nullptr;

        nc_l1a->close();

   if(l1cinput->verbose)cout<<"total number or orbital records #"<<cc-1<<endl;
   number_orecords_tot=cc-1; 
   l1cfile->norb_rec=number_orecords_tot;
    // determine ini/end time indexes for asc/desc passes
    for (size_t k = 0; k < number_orecords_tot - 1; k++) {
        if (orb_lat_tot[k + 1] > orb_lat_tot[k])
            orb_dir_tot[k] = 1;  // ascending
        else
            orb_dir_tot[k] = 0;  // descending
    }

    // find time limits for half-orbits
    int swt = 1;
    int tlimits_swt[2][2] = {{-1, -1}, {-1, -1}};  // 2 swath x time limits ini/end
    tlimits_swt[swt - 1][0] = 0;                       // only the record indexes

    for (size_t k = 0; k < number_orecords_tot - 1; k++) {
        if (orb_dir_tot[k] != orb_dir_tot[k + 1]) {
            if (swt == 1) {
                tlimits_swt[swt - 1][1] = k;  // this time is a 'local' end cause may continue circling the
                                              // earth after ascending or descending
                swt++;
           //     tlimits_swt[swt - 1][0] = k + 1;
              tlimits_swt[swt - 1][0]=0;
         //   cout<<"k ----"<<tlimits_swt[swt - 1][0]<<endl;
            } else if (swt == 2) {  // another partial swath
                tlimits_swt[swt - 1][1] = k;
            
                tlimits_swt[swt - 2][1] = number_orecords_tot - 1;
                swt++;
            }
        }
        if (k == number_orecords_tot - 2 && swt == 2)
            tlimits_swt[swt - 1][1] = number_orecords_tot - 1;
    }

    nswath = swt;

    if(l1cinput->verbose)cout << "nswath...................." << nswath << "number_orecords_tot.." << number_orecords_tot << endl;

    // split half-orbits----  group lat/lon time and velocity vector
    for (size_t sw = 0; sw < 2; sw++) {
        if (nswath == 2) {
            if (sw == 0) {
                ix1 = tlimits_swt[sw][0];
                ix2=number_orecords_tot - 1;
                if(l1cinput->verbose)cout << "nswath#2......swat#1...ix1" << ix1 << "ix2.." << ix2 << endl;
            } else if (sw == 1) {
                ix3 = tlimits_swt[sw][0];
                ix4 = number_orecords_tot - 1;
                if(l1cinput->verbose)cout << "nswath2.....swat#2...ix3" << ix3 << "ix4.." << ix4 << endl;
            }
        } else if (nswath == 3) {  //>2 we have 4 limits for swath 1
            if (sw == 0) {
                ix1 = tlimits_swt[sw][0];
                ix2 = tlimits_swt[sw + 1][0] - 1;
                ix3 = tlimits_swt[sw + 1][1] + 1;
                ix4 = tlimits_swt[sw][1];
                if(l1cinput->verbose)cout << "nswath3....swat#1...ix1" << ix1 << "ix2.." << ix2 << "ix3.." << ix3 << "ix4..."
                     << ix4 << endl;
            } else if (sw == 1) {
                ix5 = tlimits_swt[sw][0];         
                ix6 = number_orecords_tot - 1;

                if(l1cinput->verbose)cout << "nswath3....swat#2...ix5" << ix5 << "ix6.." << ix6 << endl;
            }
        } else {
            if(l1cinput->verbose)cout << "number of swaths is less than 2..!! exit..." << endl;
        
        }
    }

    float tcross1 = -999., loncross1 = -999., tcross2 = -999., loncross2 = -999.;


    l1cfile->num_gridlines = 5200;//default 4800 before
    num_gridlines = l1cfile->num_gridlines;
    //--------------------------------------------------------------------------------------
    // nswath#2
    //----case type: half orbits are consecutive
    // swath1
    //--------------------------------------------------------------------------------------
    if (nswath == 2) 
    {
        if (ix1 >= 0 && ix5 < 0) 
        {
            norbs = ix2 - ix1 + 1;               
            // compute solar zenith angle for checking day/night conditions for each swath ---
            //-------------------------------------------------------
            day_mode = sunz_swt(ix1, ix2, hour_tot, day_tot, year_tot, orb_lat_tot, orb_lon_tot);

            if(l1cinput->verbose)cout << "day_mode at nswath #2 SWATH1." << day_mode << endl;

            if (day_mode == 1) 
            {
                // orbit direction
                l1cfile->orb_dir = orb_dir_tot[ix1];

                // ini and end UTC--this is swath orbital time, not interpolated time
                tswt = (double*)calloc(norbs, sizeof(double));
                latswt = (double*)calloc(norbs, sizeof(double));
                lonswt = (double*)calloc(norbs, sizeof(double));

                status = ect_swt(l1cfile, ix1, ix2, orb_time_tot, orb_lat_tot, orb_lon_tot, orb_vel_tot,
                                 grn_vel_tot, tswt, latswt, lonswt, &tcross1, &loncross1, &mov1, &mgv1);

                if(l1cinput->verbose)cout << "nswath==2 --tcross equat crossing in (s)..swath#1..." << tcross1 << "loncross1..."
                     << loncross1 << endl;

                if (latswt != nullptr)
                    free (latswt);
                latswt = nullptr;
                if (lonswt != nullptr)
                    free (lonswt);
                lonswt = nullptr;


                if (status == 0) {
                    tmgv1 = (double*)calloc(l1cfile->num_gridlines, sizeof(double));
                    tmgvf = (double*)calloc(l1cfile->num_gridlines, sizeof(double));

                    swtime_swt2(1, l1cinput, l1cfile, norbs, tswt, tcross1, mgv1, tmgv1,orb_time_tot,number_orecords_tot);
                    // ini/end times for swath
                    num_gridlines = l1cfile->num_gridlines;
                    if(l1cinput->verbose)cout << "number of across bins L1C grid...#.." << l1cfile->nbinx << "l1cfile->n_views..."
                         << l1cfile->n_views << endl;

                    lat_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                    lon_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                    alt = allocate2d_float(num_gridlines, l1cfile->nbinx);
                   
                    orb_to_latlon(ix1,ix4,num_gridlines, l1cfile->nbinx, orb_time_tot, posr_tot,
                                  velr_tot, mgv1, tmgv1, tmgvf, lat_gd, lon_gd, alt,FirsTerrain);
                

                    if (tmgv1 != nullptr)
                        free (tmgv1);
                    tmgv1 = nullptr;

                    l1cfile->num_gridlines = l1cfile->num_gridlines - 1;
                  
                    if (l1cinput->grantype == 1) {
                        tswt_ini_sec = tfile_ini_sec + tmgvf[0];
                        tswt_end_sec = tfile_ini_sec + tmgvf[num_gridlines - 2];
                        create_time_swt(num_gridlines, tfile_ini_sec, tmgvf, tswt_ini_sec,tswt_end_sec,
                                        &tswt_ini, &tswt_ini_file,&tswt_mid, &tswt_end);

                        l1cfile->tswt_ini = tswt_ini;
                        l1cfile->tswt_mid = tswt_mid;
                        l1cfile->tswt_end = tswt_end;
                        l1cfile->tswt_ini_file = tswt_ini_file;
                        create_SOCEA2(1, l1cinput, l1cfile, lat_gd, lon_gd, alt,
                                      tmgvf);  // THIS IS SWATH PROCESSING
                    } 
                    else if (l1cinput->grantype == 0) 
                    {
                        //--------------------------------------------------------------
                        // granule processing---------------------------------------------
                        //-----------------------------------------------------------
                        deltasec = tmgvf[num_gridlines - 2] - tmgvf[0] + 1;
                        if(l1cinput->verbose)cout << "deltasec..swath." << deltasec << endl;

                        if (tswt != nullptr)
                            free (tswt);
                        tswt = nullptr;

                      //  numgran = 144 * 2;//1 day
                        numgran=12;
                        l1cfile->numgran = numgran;
                        gd_per_gran = round(num_gridlines / 10);  // 10 granules per half orbit
                        l1cfile->gd_per_gran = gd_per_gran;
                        if(l1cinput->verbose)cout << "estimated # of granules to be processed..." << numgran << "gd_per_gran..."
                             << gd_per_gran << "#gridlines.." << num_gridlines << endl;
                     
                     if(FirsTerrain)   write_L1C_granule2(1, l1cfile, l1cinput, tmgvf, lat_gd, lon_gd, alt,orb_time_tot);
                    } 
                    else 
                    {
                        cout << "ERROR selecting grantype, must be 0: granules or 1: "
                                "swath........................."
                             << endl;
                        exit(1);
                    }
                    //-----------------------------------------------------------------
                    //-----------------------------------------------------------------
                    if (lat_gd != nullptr)
                        free (lat_gd);
                    lat_gd = nullptr;
                    if (lon_gd != nullptr)
                        free (lon_gd);
                    lon_gd = nullptr;
                    if (tmgvf != nullptr)
                        free (tmgvf);
                    tmgvf = nullptr;
                    if (alt != nullptr)
                        free (alt);
                    alt = nullptr;
               
                }  // status =0

            }  // end day_mode
            else 
            {
                if(l1cinput->verbose)cout << "ERROR swath #1 day_mode = " << day_mode
                     << "nightime...continue to swath2 (nswath#2).............." << endl;
            }


        }  // end index x1 and x5 condition

    }  // end nswath=2

    //-------------HALF-ORBITS IN TWO NON-CONSECUTIVE PORTIONS ------
    // nswath =3
    // swath1
    //-------------------------------------------------------------
    c = 0;
    if (nswath == 3) {
        if (ix1 >= 0 && ix5 > 0) {
            norbs = ix2 - ix1 + 1;

            day_mode = sunz_swt(ix1, ix2, hour_tot, day_tot, year_tot, orb_lat_tot, orb_lon_tot);

            if(l1cinput->verbose)cout << "day_mode at nswath #3 SWATH1--segment#1: " << day_mode << endl;

            if (day_mode == 1) {  // daylight
                l1cfile->orb_dir = orb_dir_tot[ix1];

                tswt = (double*)calloc(norbs, sizeof(double));
                latswt = (double*)calloc(norbs, sizeof(double));
                lonswt = (double*)calloc(norbs, sizeof(double));

                status1 = ect_swt(l1cfile, ix1, ix2, orb_time_tot, orb_lat_tot, orb_lon_tot, orb_vel_tot,
                                  grn_vel_tot, tswt, latswt, lonswt, &tcross1, &loncross1, &mov1, &mgv1);
                if(l1cinput->verbose)cout << "nswath==3 --tcross equat crossing in (s)..swath#1.(segment #1).." << tcross1
                     << "loncross1..." << loncross1 << "orbit direction.0 is descending." << l1cfile->orb_dir
                     << "mov1.." << mov1 << "mgv1.." << mgv1 << endl;

                c = 0;
                if (tcross1 < 0.) {
                    if (tswt != nullptr)
                        free (tswt);
                    tswt = nullptr;
                    if (latswt != nullptr)
                        free (latswt);
                    latswt = nullptr;
                    if (lonswt != nullptr)
                        free (lonswt);
                    lonswt = nullptr;

                    norbs = ix4 - ix3 + 1;

                    tswt = (double*)calloc(norbs, sizeof(double));
                    latswt = (double*)calloc(norbs, sizeof(double));
                    lonswt = (double*)calloc(norbs, sizeof(double));

                    tcross1 = -999, loncross1 = -999.;

                    l1cfile->orb_dir = orb_dir_tot[ix3];

                    day_mode = sunz_swt(ix3, ix4, hour_tot, day_tot, year_tot, orb_lat_tot, orb_lon_tot);

                    if(l1cinput->verbose)cout << "day_mode at nswath #3 SWATH1--segment#2: " << day_mode << endl;

                    status2 = ect_swt(l1cfile, ix3, ix4, orb_time_tot, orb_lat_tot, orb_lon_tot, orb_vel_tot,
                                      grn_vel_tot, tswt, latswt, lonswt, &tcross1, &loncross1, &mov1, &mgv1);
                    if(l1cinput->verbose)cout << "nswath==3 -swath1---tcross equat crossing in (s)..swath#1.(segment #2).."
                         << tcross1 << "loncross1..." << loncross1 << "orbit direction.0 is descending."
                         << l1cfile->orb_dir << "mov1.." << mov1 << "mgv1.." << mgv1 << endl;
                }  // end segment2

             

                if (status1 == 0 || status2 == 0 && day_mode == 1) {
                    tmgv1 = (double*)calloc(l1cfile->num_gridlines, sizeof(double));

                    if(l1cinput->verbose)cout << "#gridlines..." << l1cfile->num_gridlines << "norbs.." << norbs << endl;

                    for (int i = 0; i < l1cfile->num_gridlines; i++)
                        tmgv1[i] = -1;

                    swtime_swt2_segment(
                        1, l1cinput, l1cfile, norbs, tswt, tcross1, mgv1,
                        tmgv1,orb_time_tot,number_orecords_tot);  // number of gridlines may change here!!!!!!!!! not anymore 4000 and
                                 // assymetric around the equator so bina and binb not the same!!!
                    num_gridlines = l1cfile->num_gridlines;

                    double* tmgv1_segment = (double*)calloc(l1cfile->num_gridlines, sizeof(double));
                    tmgvf = (double*)calloc(l1cfile->num_gridlines, sizeof(double));

                    int ss = 0;
                    for (int i = 0; i < num_gridlines; i++) {
                        if (tmgv1[i] >= 0) {
                            tmgv1_segment[ss] = tmgv1[i];
                            ss++;
                        }
                    }

                    if(l1cinput->verbose)cout << "#segment gridlines.ss counter.." << ss << "equal to num_gridlines..."
                         << num_gridlines << endl;

                    if (tswt != nullptr)
                        free (tswt);
                    tswt = nullptr;
                    if (latswt != nullptr)
                        free (latswt);
                    latswt = nullptr;
                    if (lonswt != nullptr)
                        free (lonswt);
                    lonswt = nullptr;

                    if(l1cinput->verbose)cout << "number of across bins L1C grid...#.." << l1cfile->nbinx << "l1cfile->n_views..."
                         << l1cfile->n_views << endl;
                    lat_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                    lon_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                    alt = allocate2d_float(num_gridlines, l1cfile->nbinx);

                    if(status1==0){
                           ix_ini=ix1;
                           ix_end=ix2;
                           }
                    else if(status2==0){
                           ix_ini=ix3;
                           ix_end=ix4;
                           }
                    orb_to_latlon(ix_ini,ix_end,num_gridlines, l1cfile->nbinx, orb_time_tot, posr_tot,
                                  velr_tot, mgv1, tmgv1_segment, tmgvf, lat_gd, lon_gd, alt,FirsTerrain);
                

                    if (tmgv1 != nullptr)
                        free (tmgv1);
                    tmgv1 = nullptr;
                    if (tmgv1_segment != nullptr)
                        free (tmgv1_segment);
                    tmgv1_segment = nullptr;

                    l1cfile->num_gridlines = l1cfile->num_gridlines - 1;

                    if (l1cinput->grantype == 1) {
                        if(l1cinput->verbose)cout << "Processing swath--->" << endl;
                        tswt_ini_sec = tfile_ini_sec + tmgvf[0];
                        tswt_end_sec = tfile_ini_sec + tmgvf[num_gridlines - 2];

                        create_time_swt(num_gridlines, tfile_ini_sec, tmgvf, tswt_ini_sec, tswt_end_sec,
                                        &tswt_ini, &tswt_ini_file, &tswt_mid,&tswt_end);

                        l1cfile->tswt_ini = tswt_ini;
                        l1cfile->tswt_mid = tswt_mid;
                        l1cfile->tswt_end = tswt_end;
                        l1cfile->tswt_ini_file = tswt_ini_file;

                        create_SOCEA2(1, l1cinput, l1cfile, lat_gd, lon_gd, alt, tmgvf);
                    } else if (l1cinput->grantype == 0) {
                        if(l1cinput->verbose)cout << "Processing granules--->" << endl;
                        numgran = 144 * 2;
                        l1cfile->numgran = numgran;
                        gd_per_gran = round(num_gridlines / 10);
                        l1cfile->gd_per_gran = gd_per_gran;

                        if(l1cinput->verbose)cout << "estimated # of granules to be processed..." << numgran << "gd_per_gran..."
                             << gd_per_gran << "#gridlines.." << num_gridlines << endl;
                    if(FirsTerrain) write_L1C_granule2(1, l1cfile, l1cinput, tmgvf, lat_gd, lon_gd, alt,orb_time_tot);
                    } else {
                        cout << "ERROR selecting grantype, must be 0: granules or 1: "
                                "swath........................."
                             << endl;
                        exit(1);
                    }

                    //-----------------------------------------------------------------
                    //-----------------------------------------------------------------
                    if (lat_gd != nullptr)
                        free (lat_gd);
                    lat_gd = nullptr;
                    if (lon_gd != nullptr)
                        free (lon_gd);
                    lon_gd = nullptr;
                    if (tmgvf != nullptr)
                        free (tmgvf);
                    tmgvf = nullptr;
                    if (alt != nullptr)
                        free (alt);
                    alt = nullptr;
                }  // end tcross (all segments)
                else {
                  if(l1cinput->verbose)
                  {  
                      cout << "ERROR swath #1 does not cross the equator..." << endl;
                      cout << "checking swath segment #2..." << endl;
                  }
                }
            }  // end day_mode
            else {
                if(l1cinput->verbose){
                cout << "day_mode==0 (nighttime)...." << day_mode << endl;
                cout << "checking swath segment #2..." << endl;
                }
            }
        }  // end indexes x1 and x5 condition

        //-----------------------------------------------------------------------------------
        // nswath =3
        // swath 2
        //----------------------------------------------------------------------------------

        norbs = ix6 - ix5 + 1;

        day_mode = sunz_swt(ix5, ix6, hour_tot, day_tot, year_tot, orb_lat_tot, orb_lon_tot);

        if(l1cinput->verbose)cout << "day_mode at nswath #3 SWATH#2: " << day_mode << endl;

        if (day_mode == 1) {
            l1cfile->orb_dir = orb_dir_tot[ix5];

            tswt2 = (double*)calloc(norbs, sizeof(double));
            latswt2 = (double*)calloc(norbs, sizeof(double));
            lonswt2 = (double*)calloc(norbs, sizeof(double));

            status = ect_swt(l1cfile, ix5, ix6, orb_time_tot, orb_lat_tot, orb_lon_tot, orb_vel_tot,
                             grn_vel_tot, tswt2, latswt2, lonswt2, &tcross2, &loncross2, &mov2, &mgv2);
            if(l1cinput->verbose)cout << "nswath==3 -swath2--tcross equat crossing in (s)..swath#2..." << tcross2 << "loncross2..."
                 << loncross2 << "orbit direction.0 is descending." << l1cfile->orb_dir << "mov2.." << mov2
                 << "mgv2.." << mgv2 << endl;

            if (latswt2 != nullptr)
                free (latswt2);
            latswt2 = nullptr;
            if (lonswt2 != nullptr)
                free (lonswt2);
            lonswt2 = nullptr;



            if (status == 0 && day_mode == 1) {//crossing and daytime
                tmgv2 = (double*)calloc(l1cfile->num_gridlines, sizeof(double));
                tmgvf2 = (double*)calloc(l1cfile->num_gridlines, sizeof(double));

                swtime_swt2(2, l1cinput, l1cfile, norbs, tswt2, tcross2, mgv2, tmgv2,orb_time_tot,number_orecords_tot);
                num_gridlines = l1cfile->num_gridlines;

                if (tswt2 != nullptr)
                    free (tswt2);
                tswt2 = nullptr;

                if(l1cinput->verbose)cout << "number of across bins L1C grid...#.." << l1cfile->nbinx << "l1cfile->n_views..."
                     << l1cfile->n_views << endl;
                lat_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                lon_gd = allocate2d_float(num_gridlines, l1cfile->nbinx);
                alt = allocate2d_float(num_gridlines, l1cfile->nbinx);

                orb_to_latlon(ix5,ix6,num_gridlines, l1cfile->nbinx, orb_time_tot, posr_tot,
                              velr_tot, mgv2, tmgv2, tmgvf2, lat_gd, lon_gd, alt,FirsTerrain);
           
               
                if (tmgv2 != nullptr)
                    free (tmgv2);
                tmgv2 = nullptr;

                l1cfile->num_gridlines = l1cfile->num_gridlines - 1;

                if (l1cinput->grantype == 1) {
                    tswt_ini_sec = tfile_ini_sec + tmgvf2[0];
                    tswt_end_sec = tfile_ini_sec + tmgvf2[num_gridlines - 2];

                    create_time_swt(num_gridlines, tfile_ini_sec, tmgvf2, tswt_ini_sec, tswt_end_sec,
                                    &tswt_ini, &tswt_ini_file, &tswt_mid,&tswt_end);

                    l1cfile->tswt_ini = tswt_ini;
                    l1cfile->tswt_mid = tswt_mid;
                    l1cfile->tswt_end = tswt_end;
                    l1cfile->tswt_ini_file = tswt_ini_file;
                    create_SOCEA2(2, l1cinput, l1cfile, lat_gd, lon_gd, alt, tmgvf2);
                } else if (l1cinput->grantype == 0) {
                    deltasec = tmgvf2[num_gridlines - 2] - tmgvf2[0] + 1;
                    if(l1cinput->verbose)cout << "deltasec..swath." << deltasec << endl;
                    numgran = 144 * 2;

                    l1cfile->numgran = numgran;
                    gd_per_gran = round(num_gridlines / 10);  // 10 granules per half orbit
                    l1cfile->gd_per_gran = gd_per_gran;

                    if(l1cinput->verbose)cout << "estimated # of granules to be processed..." << numgran << "gd_per_gran..."
                         << gd_per_gran << "#gridlines.." << num_gridlines << endl;

                 if(FirsTerrain) write_L1C_granule2(2, l1cfile, l1cinput, tmgvf2, lat_gd, lon_gd, alt,orb_time_tot);
                } else {
                    cout << "ERROR selecting grantype, must be 0: granules or 1: "
                            "swath........................."
                         << endl;
                    exit(1);
                }
                if (lat_gd != nullptr)
                    free (lat_gd);
                lat_gd = nullptr;
                if (lon_gd != nullptr)
                    free (lon_gd);
                lon_gd = nullptr;
                if (tmgvf2 != nullptr)
                    free (tmgvf2);
                tmgvf2 = nullptr;
                if (alt != nullptr)
                    free (alt);
                alt = nullptr;
            }  // status 0
            else {
              if(l1cinput->verbose)  cout << "ERROR swath #2 does not cross the equator..NO L1C grid for swath #2" << endl;
            }

        }  // end day_mode==1
        else {
            cout << "day_mode = 0 nightime...no L1C grid produced--exit (swath#2 ---nswath#3).............."
                 << day_mode << endl;
  //          exit(1);
        }
    }  // end nswath =3
 }  // end big loop files HKT batch processing
 
    free (year_tot);
    free (day_tot);
    free (hour_tot);
    free (orb_time_tot);
    free (orb_lon_tot);
    free (orb_lat_tot);
    free (orb_dir_tot);
    delete[] (posr_tot);
    delete[] (velr_tot);
    free (orb_vel_tot);
    free (grn_vel_tot);

    return 0;
}


int L1C::create_SOCEA2(int swtd, L1C_input* l1cinput, l1c_filehandle* l1cfile, float** lat_gd, float** lon_gd,
                       float** altitude, double* time_nad) {
    string tswt_ini, tswt_end, tswt_ini_file, date_created;
    int asc_mode = -1;
    std::string extstr = ".nc";
    string dirstr, prodstr;
    const char* filename_lt;
    int Nwest = -1, Neast = -1, Ngring = -1, midix = -1, dp = -1;
    int p, ix = -1;
    float latemp = -1, lontemp1 = -1, lontemp2 = -1, dlat_gd = -1, dlon_gd = -1, dlat20 = -1, dlon20 = -1,
          lon360 = -1;
    int re = -1, rw = -1;
    int32_t iyear=0,iday=0,msec=0;
    string datezerotime,yearstr,monstr,daystr;
    int16_t syear, smon, sday;
    double tdate_ini,secs;
    int twodays=0;

    if(l1cinput->verbose)cout << "Creating SOCEA proj for the whole orbit......." << endl;

    if (l1cinput->sensor == 34) {
        //"SPEXONE";
        l1cfile->nbinx = 29;
        midix = 14;
        //    NVIEWS=l1cfile->n_views;
        //    NBANDS=400;
    } else if (l1cinput->sensor == 30 || l1cinput->sensor == 31) {  // OCIS/OCIS
        //"OCI";
        l1cfile->nbinx = 519;
        midix = 259;
        //    NVIEWS=2;
        //    NBANDS=249;
    } else if (l1cinput->sensor == 35) {
        //"HARP2";
        l1cfile->nbinx = 457;
        midix = 228;
        //     NVIEWS=l1cfile->n_views;
        //     NBANDS=4;
    } else {
        cout << "sensor by default is OCI option 2....." << endl;
        //"OCI";
        l1cfile->nbinx = 519;
        midix = 259;
        //    NVIEWS=2;
        //    NBANDS=249;
    }

    for(int gd=0;gd<l1cfile->num_gridlines;gd++) {
        if(time_nad[gd]<0)
            twodays=1;
    }
    // time nadir in seconds

    tswt_ini_file = l1cfile->tswt_ini_file;
    prodstr = "PACE." + tswt_ini_file + ".L1C" + extstr;

    l1cfile->gridname = prodstr.c_str();

    filename_lt = prodstr.c_str();
    char* gridchar = strdup(filename_lt);
    string l1c_str = filename_lt;

    NcFile* nc_output;
    try {
        nc_output = new NcFile(filename_lt, NcFile::replace);
    } catch (NcException& e) {
        e.what();
        cerr << "l1cgen l1c_pflag= 5 : producing L1C grid: " + l1c_str << endl;
        exit(1);
    }

    bin_str binl1c;
    binl1c.sensor = l1cinput->sensor;
    binl1c.date_mid_grid=l1cfile->tswt_mid;
    meta_l1c_grid(gridchar, &binl1c, l1cfile->num_gridlines, nc_output);
    // gobal attrs--
    asc_mode = l1cfile->orb_dir;
    if (asc_mode == 1)
        dirstr = "Ascending";
    else if (asc_mode == 0)
        dirstr = "Descending";
    tswt_ini = l1cfile->tswt_ini;
    tswt_end = l1cfile->tswt_end;

    nc_output->putAtt("processing_version", l1cinput->pversion);
    nc_output->putAtt("history", l1cinput->history);
    nc_output->putAtt("product_name", prodstr);
    nc_output->putAtt("startDirection", dirstr);
    nc_output->putAtt("endDirection", dirstr);
    nc_output->putAtt("time_coverage_start", tswt_ini);
    nc_output->putAtt("time_coverage_end", tswt_end);

    isodate2ydmsec((char *) l1cfile->tswt_mid.c_str(), &iyear, &iday, &msec);
    double dist_es=esdist_(&iyear,&iday,&msec);
    nc_output->putAtt("sun_earth_distance", ncFloat, dist_es);
    if(l1cinput->verbose)cout << "sun_earth_distance -- mid gridline" <<dist_es<<endl;
    
    tdate_ini = isodate2unix(l1cinput->start_time); 
                  unix2ymds(tdate_ini, &syear, &smon, &sday, &secs);    
    yearstr=std::to_string(syear);
    monstr=std::to_string(smon);
    if(monstr.size() == 1)
        monstr = "0" + monstr;
    daystr=std::to_string(sday);
    if(daystr.size() == 1)
        daystr = "0" + daystr;
    datezerotime="seconds since " + yearstr + "-" + monstr + "-" + daystr;
    // vars
    NcGroup ba_grp = nc_output->getGroup("bin_attributes");
    NcVar v1 = ba_grp.getVar("nadir_view_time");

    if(twodays) { 
        for(int gd=0;gd<l1cfile->num_gridlines;gd++) {    
            time_nad[gd]+=3600*24;      
        }
    }

    v1.putVar(&time_nad[0]);
    v1.putAtt("units",datezerotime);
    //   v1=ba_grp.getVar("view_time_offsets");
    //   v1.putVar(&time_off[0][0][0]);
    NcGroup geo_grp = nc_output->getGroup("geolocation_data");
    v1 = geo_grp.getVar("latitude");
    v1.putVar(&lat_gd[0][0]);
    v1 = geo_grp.getVar("longitude");
    v1.putVar(&lon_gd[0][0]);
    v1 = geo_grp.getVar("height");
    v1.putVar(&altitude[0][0]);

    // GRING-------
    // determine the number of GCpoint indexes
    // default 6 for the swath sides + number of coordinates every 20 degrees latitude

    // ascending pass
    if (l1cfile->orb_dir == 1) {
        Nwest = round((lat_gd[l1cfile->num_gridlines - 1][0] - lat_gd[0][0]) / 20);
        Neast = round(
            (lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1] - lat_gd[0][l1cfile->nbinx - 1]) / 20);
    } else  // descending
    {
        Neast = (round(lat_gd[0][0] - lat_gd[l1cfile->num_gridlines - 1][0]) / 20);
        Nwest = round(
            (lat_gd[0][l1cfile->nbinx - 1] - lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1]) / 20);
    }

    // first NGring estimate----
    Ngring = Nwest + Neast + 6;

    if (Ngring > 0) {
        float* latarr = (float*)calloc(Ngring, sizeof(float));
        float* lonarr = (float*)calloc(Ngring, sizeof(float));     
        int* p_west = (int*)calloc(Ngring, sizeof(int));
        int* p_east = (int*)calloc(Ngring, sizeof(int));

        // corners--counterclockwise and ascending pass
        if (l1cfile->orb_dir == 1) {
            if (Ngring == 6) {
                latarr[0] = lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                latarr[1] = lat_gd[l1cfile->num_gridlines - 1][midix];
                latarr[2] = lat_gd[l1cfile->num_gridlines - 1][0];
                latarr[3] = lat_gd[0][0];
                latarr[4] = lat_gd[0][midix];
                latarr[5] = lat_gd[0][l1cfile->nbinx - 1];

                lonarr[0] = lon_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                lonarr[1] = lon_gd[l1cfile->num_gridlines - 1][midix];
                lonarr[2] = lon_gd[l1cfile->num_gridlines - 1][0];
                lonarr[3] = lon_gd[0][0];
                lonarr[4] = lon_gd[0][midix];
                lonarr[5] = lon_gd[0][l1cfile->nbinx - 1];
            } else {
                latarr[0] = lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                latarr[1] = lat_gd[l1cfile->num_gridlines - 1][midix];
                latarr[2] = lat_gd[l1cfile->num_gridlines - 1][0];

                lonarr[0] = lon_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                lonarr[1] = lon_gd[l1cfile->num_gridlines - 1][midix];
                lonarr[2] = lon_gd[l1cfile->num_gridlines - 1][0];

                latemp = latarr[2];
                latemp -= 20;
                p = 1;
                rw = 0;
                // west side
                while (latemp > lat_gd[0][0]) {
                    latarr[2 + p] = latemp;
                     if(l1cinput->verbose)cout << "p west--" << p << "point# = " << 3 + p << "lat20 = " << latemp << endl;
                    p_west[rw] = 3 + p;
                    latemp -= 20;
                    p++;
                    rw++;
                }

                p--;

                latarr[2 + p + 1] = lat_gd[0][0];
                latarr[2 + p + 2] = lat_gd[0][midix];
                latarr[2 + p + 3] = lat_gd[0][l1cfile->nbinx - 1];

                lonarr[2 + p + 1] = lon_gd[0][0];
                lonarr[2 + p + 2] = lon_gd[0][midix];
                lonarr[2 + p + 3] = lon_gd[0][l1cfile->nbinx - 1];

                latemp = latarr[2 + p + 3];
                latemp += 20;
                p++;

                int c = 1;
                re = 0;
                // east side
                while (latemp < lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1]) {
                    latarr[5 + p] = latemp;
                     if(l1cinput->verbose)cout << "p east--" << c << "point# = " << 6 + p << "lat20 = " << latemp << endl;
                    p_east[re] = 6 + p;
                    latemp += 20;
                    p++;
                    c++;
                    re++;
                }

                p--;

                 if(l1cinput->verbose){cout << "mid points west.." << rw << "mid points east.." << re << endl;
                cout << "# points in GRING = " << 6 + p << "# mid points west--" << rw
                     << "# mid points east--" << re << endl;
                 }
                if (Ngring == 6)
                    dp = 6;
                else
                    dp = rw + re + 6;
                // west
                for (int i = 0; i < rw; i++) {
                    ix = p_west[i] - 1;
                     if(l1cinput->verbose)           cout<<"lat:"<<latarr[ix]<<"pwest--"<<p_west[i]<<endl;
                    for (int row = 0; row < l1cfile->num_gridlines - 1; row++) {
                        if (latarr[ix] > lat_gd[row][0] && latarr[ix] <= lat_gd[row + 1][0]) {
                             if(l1cinput->verbose)       cout<<"mid point west# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;

                            if (lon_gd[row][0] < 0.)
                                lontemp1 = lon_gd[row][0] + 360;
                            else
                                lontemp1 = lon_gd[row][0];
                            if (lon_gd[row + 1][0] < 0.)
                                lontemp2 = lon_gd[row + 1][0] + 360;
                            else
                                lontemp2 = lon_gd[row + 1][0];

                            dlat_gd = abs(lat_gd[row + 1][0] - lat_gd[row][0]);
                            dlon_gd = abs(lontemp1 - lontemp2);
                            dlat20 = abs(latarr[ix] - lat_gd[row + 1][0]);
                            dlon20 = dlat20 * dlon_gd / dlat_gd;

                            if (lontemp1 > lontemp2)
                                lon360 = lontemp2 + dlon20;
                            else
                                lon360 = lontemp2 - dlon20;
                            if (lon360 > 180)
                                lon360 = lon360 - 360.;
                            lonarr[ix] = lon360;
                             if(l1cinput->verbose)   cout<<"lon_gd row+1.."<<lon_gd[row+1][0]<<"lon_gd row.."<<lon_gd[row][0]<<"lon cring.."<<lonarr[ix]<<endl;
                            break;
                        }
                    }
                }
                // east
                for (int i = 0; i < re; i++) {
                    ix = p_east[i] - 1;
                     if(l1cinput->verbose)cout << "lat:" << latarr[ix] << "peast--" << p_east[i] << endl;
                    for (int row = 0; row < l1cfile->num_gridlines - 1; row++) {
                        if (latarr[ix] > lat_gd[row][l1cfile->nbinx - 1] &&
                            latarr[ix] <= lat_gd[row + 1][l1cfile->nbinx - 1]) {
                             if(l1cinput->verbose)  cout<<"mid point east# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;

                            if (lon_gd[row][l1cfile->nbinx - 1] < 0.)
                                lontemp1 = lon_gd[row][l1cfile->nbinx - 1] + 360;
                            else
                                lontemp1 = lon_gd[row][l1cfile->nbinx - 1];
                            if (lon_gd[row + 1][l1cfile->nbinx - 1] < 0.)
                                lontemp2 = lon_gd[row + 1][l1cfile->nbinx - 1] + 360;
                            else
                                lontemp2 = lon_gd[row + 1][l1cfile->nbinx - 1];

                            dlat_gd =
                                abs(lat_gd[row + 1][l1cfile->nbinx - 1] - lat_gd[row][l1cfile->nbinx - 1]);
                            dlon_gd = abs(lontemp1 - lontemp2);
                            dlat20 = abs(latarr[ix] - lat_gd[row][l1cfile->nbinx - 1]);
                            dlon20 = dlat20 * dlon_gd / dlat_gd;
                            if (lontemp1 > lontemp2)
                                lon360 = lontemp1 - dlon20;
                            else
                                lon360 = lontemp1 + dlon20;
                            if (lon360 > 180)
                                lon360 = lon360 - 360.;
                            lonarr[ix] = lon360;
                             if(l1cinput->verbose) cout<<"lon_gd row+1.."<<lon_gd[row+1][l1cfile->nbinx-1]<<"lon_gd  row.."<<lon_gd[row][l1cfile->nbinx-1]<<"lon cring.."<<lonarr[ix]<<endl;

                            break;
                        }
                    }
                }
            }
        } else  // descending orb
        {
            if (Ngring == 6) {
                latarr[0] = lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                latarr[1] = lat_gd[l1cfile->num_gridlines - 1][midix];
                latarr[2] = lat_gd[l1cfile->num_gridlines - 1][0];
                latarr[3] = lat_gd[0][0];
                latarr[4] = lat_gd[0][midix];
                latarr[5] = lat_gd[0][l1cfile->nbinx - 1];

                lonarr[0] = lon_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                lonarr[1] = lon_gd[l1cfile->num_gridlines - 1][midix];
                lonarr[2] = lon_gd[l1cfile->num_gridlines - 1][0];
                lonarr[3] = lon_gd[0][0];
                lonarr[4] = lon_gd[0][midix];
                lonarr[5] = lon_gd[0][l1cfile->nbinx - 1];
            } else {
                latarr[0] = lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                latarr[1] = lat_gd[l1cfile->num_gridlines - 1][midix];
                latarr[2] = lat_gd[l1cfile->num_gridlines - 1][0];

                lonarr[0] = lon_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1];
                lonarr[1] = lon_gd[l1cfile->num_gridlines - 1][midix];
                lonarr[2] = lon_gd[l1cfile->num_gridlines - 1][0];

                latemp = latarr[2];
                latemp += 20;
                p = 1;
                int rw = 0;
                // west side
                while (latemp < lat_gd[0][0]) {
                    latarr[2 + p] = latemp;
                     if(l1cinput->verbose)cout << "p west--" << p << "point# = " << 3 + p << "lat20 = " << latemp << endl;
                    p_west[rw] = 3 + p;
                    latemp += 20;
                    p++;
                    rw++;
                }

                p--;

                latarr[2 + p + 1] = lat_gd[0][0];
                latarr[2 + p + 2] = lat_gd[0][midix];
                latarr[2 + p + 3] = lat_gd[0][l1cfile->nbinx - 1];

                lonarr[2 + p + 1] = lon_gd[0][0];
                lonarr[2 + p + 2] = lon_gd[0][midix];
                lonarr[2 + p + 3] = lon_gd[0][l1cfile->nbinx - 1];
                latemp = latarr[2 + p + 3];
                latemp -= 20;
                p++;

                int c = 1;
                int re = 0;
                // east side
                while (latemp > lat_gd[l1cfile->num_gridlines - 1][l1cfile->nbinx - 1]) {
                    latarr[5 + p] = latemp;
                     if(l1cinput->verbose)cout << "p east--" << c << "point# = " << 6 + p << "lat20 = " << latemp << endl;
                    p_east[re] = 6 + p;
                    latemp -= 20;
                    p++;
                    c++;
                    re++;
                }

                p--;

                 if(l1cinput->verbose){
                cout << "mid points west.." << rw << "mid points east.." << re << endl;
                cout << "# points in GRING = " << 6 + p << "# mid points west--" << rw                    
                     << "# mid points east--" << re << endl;
                 }
                if (Ngring == 6)
                    dp = 6;
                else
                    dp = rw + re + 6;

                // west side
                for (int i = 0; i < rw; i++) {
                    ix = p_west[i] - 1;
                     if(l1cinput->verbose)      cout<<"lat:"<<latarr[ix]<<"pwest--"<<p_west[i]<<endl;
                    for (int row = 0; row < l1cfile->num_gridlines - 1; row++) {
                        if (latarr[ix] <= lat_gd[row][0] && latarr[ix] > lat_gd[row + 1][0]) {
                             if(l1cinput->verbose)   cout<<"mid point west# = "<<i+1<<"found index between #row= "<<row+1<<"and   row ="<<row+2<<endl;
                            if (lon_gd[row][0] < 0.)
                                lontemp1 = lon_gd[row][0] + 360;
                            else
                                lontemp1 = lon_gd[row][0];
                            if (lon_gd[row + 1][0] < 0.)
                                lontemp2 = lon_gd[row + 1][0] + 360;
                            else
                                lontemp2 = lon_gd[row + 1][0];

                            dlat_gd = abs(lat_gd[row][0] - lat_gd[row + 1][0]);
                            dlon_gd = abs(lontemp1 - lontemp2);
                            dlat20 = abs(latarr[ix] - lat_gd[row + 1][0]);
                            dlon20 = dlat20 * dlon_gd / dlat_gd;

                            if (lontemp1 > lontemp2)
                                lon360 = lontemp2 + dlon20;
                            else
                                lon360 = lontemp2 - dlon20;
                            if (lon360 > 180)
                                lon360 = lon360 - 360.;
                            lonarr[ix] = lon360;
                             if(l1cinput->verbose)   cout<<"lon_gd row+1.."<<lon_gd[row+1][0]<<"lon_gd row.."<<lon_gd[row][0]<<"lon cring.."<<lonarr[ix]<<endl;
                            break;
                        }
                    }
                }
                // east side
                for (int i = 0; i < re; i++) {
                    ix = p_east[i] - 1;
                     if(l1cinput->verbose)    cout<<"lat:"<<latarr[ix]<<"peast--"<<p_east[i]<<endl;
                    for (int row = 0; row < l1cfile->num_gridlines - 1; row++) {
                        if (latarr[ix] <= lat_gd[row][l1cfile->nbinx - 1] &&
                            latarr[ix] > lat_gd[row + 1][l1cfile->nbinx - 1]) {
                             if(l1cinput->verbose)cout<<"mid point east# = "<<i+1<<"found index between #row= "<<row+1<<"and row ="<<row+2<<endl;
                            if (lon_gd[row][l1cfile->nbinx - 1] < 0.)
                                lontemp1 = lon_gd[row][l1cfile->nbinx - 1] + 360;
                            else
                                lontemp1 = lon_gd[row][l1cfile->nbinx - 1];
                            if (lon_gd[row + 1][l1cfile->nbinx - 1] < 0.)
                                lontemp2 = lon_gd[row + 1][l1cfile->nbinx - 1] + 360;
                            else
                                lontemp2 = lon_gd[row + 1][l1cfile->nbinx - 1];

                            dlat_gd =
                                abs(lat_gd[row][l1cfile->nbinx - 1] - lat_gd[row + 1][l1cfile->nbinx - 1]);
                            dlon_gd = abs(lontemp1 - lontemp2);
                            dlat20 = abs(latarr[ix] - lat_gd[row][l1cfile->nbinx - 1]);
                            dlon20 = dlat20 * dlon_gd / dlat_gd;
                            if (lontemp1 > lontemp2)
                                lon360 = lontemp1 - dlon20;
                            else
                                lon360 = lontemp1 + dlon20;
                            if (lon360 > 180)
                                lon360 = lon360 - 360.;
                            lonarr[ix] = lon360;
                             if(l1cinput->verbose)  cout<<"lon_gd row+1.."<<lon_gd[row+1][l1cfile->nbinx-1]<<"lon_gd row.."<<lon_gd[row][l1cfile->nbinx-1]<<"lon cring.."<<lonarr[ix]<<endl;

                            break;
                        }
                    }
                }

            }  // else GRING>6

        }  // end descending

        // sequence- GRING

       string onelat,onelon,onecoor,firstcoor;
                //GRing params
       string gring_polygon,gring_crs,gring_latmin,gring_latmax,gring_lonmin,gring_lonmax;

       for (int s = 0; s < dp; s++) {                 
                    onelat=to_string(latarr[s]);
                    onelon=to_string(lonarr[s]);
                    if(s==0)
                    {
                        onecoor="POLYGON(("+onelat+" "+onelon+",";                        
                        firstcoor=onelat+" "+onelon;
                    }
                    else onecoor=onelat+" "+onelon+","; 
                    gring_polygon+=onecoor;
                }
                //last coor
        gring_polygon+=firstcoor+"))";
        float latpmin=999,latpmax=-999,lonpmin=999,lonpmax=-999;
        for(int row=0;row<l1cfile->num_gridlines;row++)
                {
                for(int col=0;col<l1cfile->nbinx;col++)
                    {
                      if(lat_gd[row][col]<latpmin &&  lat_gd[row][col]!=BAD_FLT) latpmin=lat_gd[row][col];
                      if(lat_gd[row][col]>latpmax &&  lat_gd[row][col]!=BAD_FLT) latpmax=lat_gd[row][col];
                      if(lon_gd[row][col]<lonpmin &&  lon_gd[row][col]!=BAD_FLT) lonpmin=lon_gd[row][col];
                      if(lon_gd[row][col]>lonpmax &&  lon_gd[row][col]!=BAD_FLT) lonpmax=lon_gd[row][col];
                    }}

        gring_crs="EPSG:4326";
        nc_output->putAtt("geospatial_bounds", gring_polygon);
        nc_output->putAtt("geospatial_bounds_crs", gring_crs);

                //algo to find min/max coordinates inside the gring polygon
        nc_output->putAtt("geospatial_lat_min", ncFloat, latpmin);
        nc_output->putAtt("geospatial_lat_max", ncFloat, latpmax);
        nc_output->putAtt("geospatial_lon_min", ncFloat, lonpmin);
        nc_output->putAtt("geospatial_lon_max", ncFloat, lonpmax);

        free (latarr);
        free (lonarr);
        free (p_west);
        free (p_east);

    } else {
        cout << "ERROR EXTRACTING GRING coordinates!!-----" << endl;
        exit(1);
    }

    nc_output->close();
    return 0;
}

// fred/me implementation OK 3/22/2023
int L1C::search_l1cgen(L1C_input* l1cinput, l1c_str* l1cstr, l1c_filehandle* l1cfile, short** gdindex) {
    int32_t num_gridlines, nbinx;
    float gnvm;
    float gnvec[3], gvec[3], bvec[3];
    int flag_out = -1;
    size_t pix;
    int32_t i;
    size_t num_pixels;
    int irow = -1, col = -1;
    float bmcm, bm = 100;
    double db;
    float c1, c2, c3;
    double fudge = 0.00001, dotprod, dot_firstline, dot_lastline;
    flag_out = -1;

    num_gridlines = l1cfile->num_gridlines;
    nbinx = l1cfile->nbinx;

    num_pixels = l1cfile->npix;

    db = (l1cinput->grid_resolution) / 6371 / 2;  // Half of bin size in radians

    // big loop
    // dot product
    for (pix = 0; pix < num_pixels; pix++) {
        bvec[0] = cos(l1cstr->lonpix[pix] * M_PI / 180) * cos(l1cstr->latpix[pix] * M_PI / 180);
        bvec[1] = sin(l1cstr->lonpix[pix] * M_PI / 180) * cos(l1cstr->latpix[pix] * M_PI / 180);
        bvec[2] = sin(l1cstr->latpix[pix] * M_PI / 180);

        for (i = 0; i < num_gridlines; i++) {
            if (l1cstr->latpix[pix] > 90 || l1cstr->latpix[pix] < -90 || l1cstr->lonpix[pix] < -180 ||
                l1cstr->lonpix[pix] > 180) {
             /*   cout << "ERROR in search_l1cgen --latitude longitude pixel out of the boundaries.."
                     << "latpix>90 or <-90.." << l1cstr->latpix[pix] << "lonpix>180 or <-180.."
                     << l1cstr->lonpix[pix] << endl;*/
                return 110;
            }
            // normal vectors for L1C rows
            if (l1cfile->lat_gd[i][nbinx - 1] > 90 || l1cfile->lat_gd[i][nbinx - 1] < -90 ||
                l1cfile->lon_gd[i][nbinx - 1] < -180 || l1cfile->lon_gd[i][nbinx - 1] > 180) {
                cout << "lat lon out of the boundaries.." << endl;
                exit(1);
            }
            if (l1cfile->lat_gd[i][0] > 90 || l1cfile->lat_gd[i][0] < -90 || l1cfile->lon_gd[i][0] < -180 ||
                l1cfile->lon_gd[i][0] > 180) {
                cout << "lat lon out of the boundaries.." << endl;
                exit(1);
            }

            gnvec[0] = sin(l1cfile->lon_gd[i][nbinx - 1] * M_PI / 180) *
                           cos(l1cfile->lat_gd[i][nbinx - 1] * M_PI / 180) *
                           sin(l1cfile->lat_gd[i][0] * M_PI / 180) -
                       sin(l1cfile->lat_gd[i][nbinx - 1] * M_PI / 180) *
                           sin(l1cfile->lon_gd[i][0] * M_PI / 180) * cos(l1cfile->lat_gd[i][0] * M_PI / 180);
            gnvec[1] = sin(l1cfile->lat_gd[i][nbinx - 1] * M_PI / 180) *
                           cos(l1cfile->lon_gd[i][0] * M_PI / 180) * cos(l1cfile->lat_gd[i][0] * M_PI / 180) -
                       cos(l1cfile->lon_gd[i][nbinx - 1] * M_PI / 180) *
                           cos(l1cfile->lat_gd[i][nbinx - 1] * M_PI / 180) *
                           sin(l1cfile->lat_gd[i][0] * M_PI / 180);
            gnvec[2] = cos(l1cfile->lon_gd[i][nbinx - 1] * M_PI / 180) *
                           cos(l1cfile->lat_gd[i][nbinx - 1] * M_PI / 180) *
                           sin(l1cfile->lon_gd[i][0] * M_PI / 180) * cos(l1cfile->lat_gd[i][0] * M_PI / 180) -
                       sin(l1cfile->lon_gd[i][nbinx - 1] * M_PI / 180) *
                           cos(l1cfile->lat_gd[i][nbinx - 1] * M_PI / 180) *
                           cos(l1cfile->lon_gd[i][0] * M_PI / 180) * cos(l1cfile->lat_gd[i][0] * M_PI / 180);

            // vector norm
            gnvm = sqrt(gnvec[0] * gnvec[0] + gnvec[1] * gnvec[1] + gnvec[2] * gnvec[2]);
            if (isnan(gnvm) == 1) {
                cout << "NAN value for gnvm.." << endl;
                exit(1);
            }
            if (gnvm == 0) {
                cout << "ERROR gnvm == 0--- WE CANT NORMALIZE..." << endl;
                exit(1);
            }
            // normalization
            gnvec[0] = gnvec[0] / gnvm;
            gnvec[1] = gnvec[1] / gnvm;
            gnvec[2] = gnvec[2] / gnvm;

            // for each pixels
            // dot prod, orbital normaliz and transposed by pixel vector
            dotprod = gnvec[0] * bvec[0] + gnvec[1] * bvec[1] + gnvec[2] * bvec[2];

            if (i == 0) {
                dot_firstline = dotprod;
            }
            if (i == num_gridlines - 1) {
                dot_lastline = dotprod;
            }

            if (dotprod - fudge <= db && dotprod + fudge > -db) {
                gdindex[pix][0] = i + 1;  // first found
            }

        }  // end lines
           // for each pixels
        if (dot_firstline <= db && dot_lastline > -db) {

            // find column
            irow = gdindex[pix][0] - 1;
            if (irow < 0) {
                cout << "ERROR icol in search_l1c..."
                     << "at pix#.." << pix + 1 << "and irow#.." << irow << endl;
                exit(1);
            }

            for (int j = 0; j < nbinx; j++) {
                gvec[0] =
                    cos(l1cfile->lon_gd[irow][j] * M_PI / 180) * cos(l1cfile->lat_gd[irow][j] * M_PI / 180);
                gvec[1] =
                    sin(l1cfile->lon_gd[irow][j] * M_PI / 180) * cos(l1cfile->lat_gd[irow][j] * M_PI / 180);
                gvec[2] = sin(l1cfile->lat_gd[irow][j] * M_PI / 180);

                c1 = bvec[0] - gvec[0];
                c2 = bvec[1] - gvec[1];
                c3 = bvec[2] - gvec[2];

                bmcm = sqrt(c1 * c1 + c2 * c2 + c3 * c3);  ////bmcm only checked one pixel!!I
                if (bmcm < bm) {
                    bm = bmcm;
                    col = j + 1;
                }
            }
            if (col < 1) {
                cout << "ERROR col in search_l1c..."
                     << "at pix#.." << pix + 1 << "and row#.." << irow + 1 << endl;
                exit(1);
            }

            gdindex[pix][1] = col;
            bm = 100;
            col = -1;
        } else {
            gdindex[pix][0] = -1;  // actually these are row/col not irow/icol like in my searching algos!!!
            gdindex[pix][1] = -1;
        }

    }  // end pixels

    flag_out = -1;
    for (pix = 0; pix < num_pixels; pix++) {
        if (gdindex[pix][0] > 0 && gdindex[pix][1] > 0) {
            flag_out = 0;
        } else {
             if(l1cinput->verbose)    cout<<"THIS LINE WILL BE SKIPPED -- NO PIXELS BINNED.............."<<endl;
        }
    }

    if (flag_out == 0)
        return 0;
    else
        return 1;
}


// version after Fred proj
int32_t L1C::openL1Cgrid3(l1c_str* l1cstr, l1c_filehandle* l1cfile, L1C_input* l1cinput) {
    int32_t num_gridlines, nbinx;  // number of gridlines to be processed
    std::string str;
    const char* ptstr;
    int32_t NY = -1, NX = -1;
    std::string ifile_str;
    string gridname, azeast_name;
    string tswt_ini_file;
    std::string fname_out, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, extstr, granstr,
        timestr, azstr, missionstr, ofilestr;
    NcFile* nc_l1c;

    // char tmp[256];
    //        getcwd(tmp, 256);
    //        pathstr=tmp;

    if (l1cinput->projection == 0) {  // socea
        gridname = l1cinput->l1c_grid;
         if(l1cinput->verbose)cout << "projection type is SOCEA!!!........for grid file..:" << gridname << endl;
    }
    if (l1cinput->projection == 1) {  // fixed bearing projection ----
        gridname = "/accounts/mamontes/images/OCIS/sean/out/FB_L1Cgrid.nc";
    }

    if(l1cinput->verbose)cout << "gridname.." << gridname << endl;
    ptstr = gridname.c_str();

     if(l1cinput->verbose)cout << "Opening L1C grid ." << gridname << "for parallax-cloud corrections......." << endl;
    try {
        nc_l1c = new NcFile(ptstr, NcFile::read);
    } catch (NcException& e) {
        e.what();
        cerr << "l1cgen l1c_pflag= 8:: Failure write FULL L1C grid: " + gridname << endl;
        exit(1);
    }

    // global attributs
    NcGroupAtt i1 = nc_l1c->getAtt("instrument");
    i1.getValues(l1cfile->instrument);
    i1 = nc_l1c->getAtt("startDirection");  // NcGroupAtt is a global attr!!
    i1.getValues(l1cfile->start_dir);
    i1 = nc_l1c->getAtt("endDirection");
    i1.getValues(l1cfile->end_dir);
    i1 = nc_l1c->getAtt("time_coverage_start");
    i1.getValues(l1cfile->tswt_ini);
    i1 = nc_l1c->getAtt("time_coverage_end");
    i1.getValues(l1cfile->tswt_end);
    i1 = nc_l1c->getAtt("nadir_bin");
    i1.getValues(l1cfile->binstr);

    // dims
    NcDim ydim = nc_l1c->getDim("bins_along_track");
    NcDim xdim = nc_l1c->getDim("bins_across_track");

    //************************************************************
    // mem allocation for geolocation pointers---
    num_gridlines = ydim.getSize();
    nbinx = xdim.getSize();

    if(l1cinput->verbose)cout << "num_gridlines..." << num_gridlines << "nbinx....." << nbinx << endl;

    l1cfile->lat_gd = allocate2d_float(num_gridlines, nbinx);
    l1cfile->lon_gd = allocate2d_float(num_gridlines, nbinx);
    l1cfile->alt_gd = allocate2d_float(num_gridlines, nbinx);
    l1cfile->index_xy = allocate2d_short(num_gridlines, nbinx);
    l1cfile->lat_asort = allocate2d_float(num_gridlines, nbinx);

    NcGroup geo_grp = nc_l1c->getGroup("geolocation_data");
    NcVar v1 = geo_grp.getVar("latitude");
    v1.getVar(&l1cfile->lat_gd[0][0]);
    v1 = geo_grp.getVar("longitude");
    v1.getVar(&l1cfile->lon_gd[0][0]);
    v1 = geo_grp.getVar("altitude");
    v1.getVar(&l1cfile->alt_gd[0][0]);

    nc_l1c->close();

    NY = num_gridlines;
    NX = nbinx;
    l1cfile->num_gridlines = NY;
    l1cfile->nbinx = NX;

    //*********** SORTING LAT/LON ASCENDING AND TRACK INDEXES ****************************
    //***********************************************************************************
    vector<pair<float, int>> vp;

    // fill vector with lat_gd column--
     if(l1cinput->verbose)cout << "*********** SORTING LAT/LON GRIDPOINTS AND TRACKING INDEXES ********************" << endl;

    for (int col = 0; col < NX; col++) {
        for (int row = 0; row < NY; row++) {
            vp.push_back(make_pair(l1cfile->lat_gd[row][col] + 90., row));
        }
        // sort
        stable_sort(vp.begin(), vp.end());

        for (unsigned int i = 0; i < vp.size(); i++) {
            l1cfile->lat_asort[i][col] = vp[i].first;
            l1cfile->index_xy[i][col] = vp[i].second;
        }
        vp.clear();
    }

    return 0;
}



int32_t L1C::write_L1C_granule2(int swtd, l1c_filehandle* l1cfile, L1C_input* l1cinput, double* tmgv,
                                float** lat_gd, float** lon_gd, float** alt_gd,double *orb_time_tot) {
    int32_t num_gridlines, nbinx, NY1 = -1, NY2 = -1;
    int32_t NX, NY;
    const char* filename_lt;
    float **lat_out = nullptr, **lon_out = nullptr, **alt_out = nullptr;
    double* time_nad_out = nullptr;
    int asc_mode = l1cfile->orb_dir;
    string senstr;
    double tg_ini, tg_end,tg_mid;
    int numgran;
    double tgridline=-1, tfile_ini_sec=-1, tfile_end_sec=-1;
    int16_t* granid = nullptr;
    int32_t gransize;
    short** gdindex = nullptr;
    double** gdtime = nullptr;
    int16_t y_ini, mo_ini, d_ini, h_ini, mi_ini, y_end, mo_end, d_end, h_end, mi_end, y_mid, mo_mid, d_mid, h_mid, mi_mid,y_zero, mo_zero, d_zero,h_zero,mi_zero;
    double sec_ini, sec_end,sec_mid,sec_zero;
    string tswt_ini, tswt_end, tswt_ini_file;

    int16_t syear, smon, sday, syear2, smon2, sday2;
    double secs, secs2;
    double tgran_ini, tgran_ini_sec, tgran_end, tgran_end_sec;
    string gfull;
    int16_t gtime;
    int Ngring = -1,  dp = -1;
    int32_t iyear=0,iday=0,msec=0;    
    int twodays=0,nadir_bin_index; 
    Geodesic geod(Constants::WGS84_a(), Constants::WGS84_f());
    GeodesicLine line;
    size_t norb_rec=l1cfile->norb_rec;
    double tfile_ini_offset=-1;
  
    if (l1cinput->sensor == 34) {
        senstr = "SPEXONE";
        l1cfile->nbinx = 29;
        nadir_bin_index=14;
    
    } else if (l1cinput->sensor == 30 || l1cinput->sensor == 31) {
        senstr = "OCI";  // OCIS
        l1cfile->nbinx = 519;
        nadir_bin_index=259;
   
    } else if (l1cinput->sensor == 35) {
        senstr = "HARP2";
        l1cfile->nbinx = 457;
        nadir_bin_index=228;
    } else {
        if(l1cinput->verbose)cout << "sensor by default is OCI option 2....." << endl;
        senstr = "OCI";
        l1cfile->nbinx = 519;
        nadir_bin_index=259;
    }


    gransize = l1cfile->gransize;

    num_gridlines = l1cfile->num_gridlines;
    if(l1cinput->verbose)cout << "swath#....................................................................." << swtd
         << "asc_mode..." << asc_mode << endl;

    granid = (int16_t*)calloc(num_gridlines, sizeof(int16_t));

    numgran = l1cfile->numgran;//1 day of granules or 288x5 minutes/60 = 24 h

    nbinx = l1cfile->nbinx;

    gdtime = allocate2d_double(numgran, 2);
    gdindex = allocate2d_short(numgran, 2);

    tfile_ini_sec = l1cfile->tfile_ini_sec;
    tfile_ini_offset=tfile_ini_sec;
    tfile_end_sec = l1cfile->tfile_end_sec;
    double time_zero=tfile_ini_sec-24*3600;
    tg_ini = tfile_ini_sec;           // always unix time
    tg_end = tg_ini + gransize * 60;  // always in seconds

    // first granule----
    if (l1cinput->start_time[0] != '\0' && l1cinput->end_time[0] != '\0' && l1cinput->grantype == 0) {
        if(l1cinput->verbose)cout << "Processing L1C granules between ........................" << l1cinput->start_time
             << "and ........................." << l1cinput->end_time << endl;
        tgran_ini = isodate2unix(l1cinput->start_time);
        tgran_end = isodate2unix(l1cinput->end_time);
        unix2ymds(tgran_ini, &syear, &smon, &sday, &secs);
        unix2ymds(tgran_end, &syear2, &smon2, &sday2, &secs2);
        tgran_ini_sec = tgran_ini;

        tg_ini = tgran_ini_sec;
        tg_end = tg_ini + gransize * 60;
        tgran_end_sec = tgran_end;

        tfile_end_sec = tgran_end_sec + gransize * 60;

        unix2ymdhms(tgran_ini_sec, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
        if(l1cinput->verbose)cout << "tgran_ini_sec.."
             << "YEAR.." << y_ini << "MONTH.." << mo_ini << "DAY.." << d_ini << "HOUR.." << h_ini << "MIN.."
             << mi_ini << "SEC.." << sec_ini << endl;
        unix2ymdhms(tgran_end_sec, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);

        if (tgran_end_sec > tfile_end_sec) {
            cout << "ERROR--tgran_end_sec>tfile_end_sec...wrong command line (gran_end_time).gran_end_time "
                    "beyond swath time limits.."
                 << endl;
            exit(1);
        }

        if(l1cinput->verbose)cout << "tgran_end_sec.."
             << "YEAR.." << y_ini << "MONTH.." << mo_ini << "DAY.." << d_ini << "HOUR.." << h_ini << "MIN.."
             << mi_ini << "SEC.." << sec_ini << endl;
        unix2ymdhms(tg_ini, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
        if(l1cinput->verbose)cout << "tg_ini.."
             << "YEAR.." << y_ini << "MONTH.." << mo_ini << "DAY.." << d_ini << "HOUR.." << h_ini << "MIN.."
             << mi_ini << "SEC.." << sec_ini << endl;
        unix2ymdhms(tg_end, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
        if(l1cinput->verbose)cout << "tg_end.."
             << "YEAR.." << y_ini << "MONTH.." << mo_ini << "DAY.." << d_ini << "HOUR.." << h_ini << "MIN.."
             << mi_ini << "SEC.." << sec_ini << endl;
    } else {
         if(l1cinput->verbose)cout << "ERROR Processing L1C granules between initial and final time --"
             << "start_time.." << l1cinput->start_time << "end_time.." << l1cinput->end_time << "grantype..."
             << l1cinput->grantype << endl;
    }

    for (int i = 0; i < num_gridlines; i++) {
        granid[i] = -1;
    }

    for (int i = 0; i < numgran; i++) {
        for (size_t j = 0; j < 2; j++) {
            gdtime[i][j] = -1;
            gdindex[i][j] = -1;
        }
    }

    int neg = 0, gmin = -1, c = 0;
    if(tmgv[0]<0) twodays=1;//the granule can start later even tough the grid contains gridlines from previous day


    for (int gran = 0; gran < numgran; gran++) {
        for (int i = 0; i < num_gridlines; i++) {
            tgridline =
                tfile_ini_sec + tmgv[i];  // seconds of the day are added to second since unix time reference     
            if (tg_end > tgran_end_sec)
                tg_end = tgran_end_sec;

            // unix2ymdhms(tg_ini, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
            // cout<<"gran# "<<gran+1<<"day ini "<<d_ini<<"h ini "<<h_ini<<"mi ini "<<mi_ini<<"sec ini "<<sec_ini<<"tgini "<<tg_ini<<endl;
            //  unix2ymdhms(tg_end, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
            // cout<<"gran# "<<gran+1<<"day ini "<<d_ini<<"h ini "<<h_ini<<"mi ini "<<mi_ini<<"sec ini "<<sec_ini<<"tgini "<<tg_ini<<endl;
            // unix2ymdhms(tgridline, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
            // cout<<"gridline# "<<i+1<<"gran# "<<gran+1<<"day ini "<<d_ini<<"h ini "<<h_ini<<"mi ini "<<mi_ini<<"sec ini "<<sec_ini<<"tgini "<<tg_ini<<endl;

            if (tgridline >= tg_ini && tgridline < tg_end) {  

                if (l1cinput->start_time[0] != '\0' || l1cinput->end_time[0] != '\0') {
                    if (tg_ini >= tgran_ini_sec && tg_end <= tgran_end_sec + gransize * 60) {
                        if (gmin < 0) {  // first time of the granule or th_ini
                            gdtime[c][0] = tg_ini;
                            if (i == num_gridlines - 1)
                                gdtime[c][1] = tgridline;
                            else
                                gdtime[c][1] = tg_end;
                            gdindex[c][0] = i;
                            gmin = 1;
                        }
                        gdindex[c][1] = i;
                        if (i == num_gridlines - 1)
                            gdtime[c][1] = tgridline;
                        else
                            gdtime[c][1] = tg_end;
                    }  // gran time
                }      // command line
                else {
                    if (gmin < 0) {  
                        gdtime[c][0] = tg_ini;
                        if (i == num_gridlines - 1)
                            gdtime[c][1] = tgridline;
                        else
                            gdtime[c][1] = tg_end;
                        gdindex[c][0] = i;
                        gmin = 1;
                    }
                    gdindex[c][1] = i;
                    if (i == num_gridlines - 1)
                        gdtime[c][1] = tgridline;
                    else
                        gdtime[c][1] = tg_end;
                }
            }
      
      
            if (tmgv[i] < 0 && gran == 0) {  // previous day data
                neg++;
            }
        }  // gridlines
           // new granule----
        tg_ini = tg_end;
        tg_end = tg_ini + gransize * 60;
        gmin = -1;
        c++;
       
        if (tg_ini > tgridline) {
             if(l1cinput->verbose)cout << "gridlines with negative time..." << neg << endl;
            break;
        }
        // granule selection constraints----------------
        if (l1cinput->start_time[0] != '\0' || l1cinput->end_time[0] != '\0') {
            if (tg_ini >= tg_end || tg_ini >= tgran_end_sec || tg_end > (tgran_end_sec + gransize * 60)) {
               if(l1cfile->verbose)cout<<"ERROR : GRAN # "<<c+1<<"tg_ini >= tg_end || tg_ini >= tgran_end_sec || tg_end > (tgran_end_sec + gransize * 60"<<endl;
                break;
            }
        }

    }  // granules

    std::string timestr, missionstr, fname_out, fname_out_nopath, pathstr, monstr, daystr, yearstr, secstr,
        mistr, hstr, prodstr, gdstr, swtstr, extstr, ofilestr, dirstr1, dirstr2,datetimestr1, datetimestr2,datetimestr3,
        fdatetimestr1, date_created,datezerotime;
    std::string y_create, m_create, d_create, t_create;

    pathstr = "";
    missionstr = "PACE";
    extstr = ".nc";
    // write filenames to list----
    if(l1cinput->outlist[0]=='\0')
       strcpy(l1cinput->outlist, "l1c.tmp");
    string outxt(l1cinput->outlist);

    outxt = pathstr + outxt;

    std::ofstream outf;

    outf.open(outxt, std::ofstream::out);

    if (outf) {
         if(l1cinput->verbose)cout << "writing L1C granules to outfile..." << outxt << endl;
    } else {
        std::cerr << "output file.." << outxt << " could not be opened for writing!\n";
        return 1;
    }

    // Current date/time based on current system
    date_created = unix2isodate(now(), 'G');

    int16_t ngridlines, totlines = 0;
    // write files with granid>0 --------------------- 
    for (int gran = 0; gran < numgran; gran++) {
  

        NY1 = gdindex[gran][0];
        NY2 = gdindex[gran][1];
 




        if(twodays)
        { 
        tfile_ini_offset=time_zero;
        if(l1cinput->verbose==1) cout<<"twodays flag is ON!! : "<<twodays<<"tzero offset = "<< tfile_ini_offset<<endl;
        }   
        
                 
        if(l1cinput->verbose){
             tg_ini = gdtime[gran][0];
             tg_end = gdtime[gran][1];
             unix2ymdhms(tg_ini, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
             unix2ymdhms(tg_end, &y_end, &mo_end, &d_end, &h_end, &mi_end, &sec_end);
             cout<<"--------------------------"<<endl;
             cout<<"gran# "<<gran+1<<"day ini "<<d_ini<<"h ini "<<h_ini<<"mi ini "<<mi_ini<<"sec ini "<<sec_ini<<"tgini "<<tg_ini<<endl;
             cout<<"gran# "<<gran+1<<"day end "<<d_end<<"h end "<<h_end<<"mi end "<<mi_end<<"sec end "<<sec_end<<"tend "<<tg_end<<endl;
             tg_ini = tfile_ini_offset+orb_time_tot[0];
             tg_end = tfile_ini_sec+orb_time_tot[norb_rec-1];
             unix2ymdhms(tg_ini, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
             unix2ymdhms(tg_end, &y_end, &mo_end, &d_end, &h_end, &mi_end, &sec_end);
             cout<<"orbital time limits ----------------- "<<endl;
             cout<<"gran# "<<gran+1<<"day ini "<<d_ini<<"h ini "<<h_ini<<"mi ini "<<mi_ini<<"sec ini "<<sec_ini<<"tgini "<<tg_ini<<endl;
             cout<<"gran# "<<gran+1<<"day end "<<d_end<<"h end "<<h_end<<"mi end "<<mi_end<<"sec end "<<sec_end<<"tend "<<tg_end<<endl;
             cout<<"---------------------------"<<endl;
             }

         double torb_off=0;//5 more minutes after orb to avoid time truncation by granule
         //some orbital records from previous day but not in the crossing granule
         if(tfile_ini_offset+orb_time_tot[0]>tfile_ini_sec+orb_time_tot[norb_rec-1]) orb_time_tot[0]-=24*3600;
         if(gdtime[gran][1]>tfile_ini_sec+orb_time_tot[norb_rec-1]) torb_off=gdtime[gran][1]-tfile_ini_sec-orb_time_tot[norb_rec-1];

         if(torb_off>2*60) torb_off=0;

         if (gdtime[gran][0]>=(tfile_ini_offset+orb_time_tot[0]) && gdtime[gran][1]<=(tfile_ini_sec+orb_time_tot[norb_rec-1]+torb_off)) 
         {
            if(l1cinput->verbose)cout << "gran #..." << gran + 1 <<endl;

            tg_ini = gdtime[gran][0];
            tg_end = gdtime[gran][1];      
            tg_mid=(tg_ini+tg_end)/2;
           
            unix2ymdhms(tg_ini, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
            unix2ymdhms(tg_end, &y_end, &mo_end, &d_end, &h_end, &mi_end, &sec_end);
            unix2ymdhms(tg_mid, &y_mid, &mo_mid, &d_mid, &h_mid, &mi_mid, &sec_mid);

            if (mi_end * 60 == 0) 
            {
                gtime = ((60 * 60 + round(sec_end)) - mi_ini * 60 - round(sec_ini)) / 60;
            } 
            else 
            {
                gtime = ((mi_end * 60 + round(sec_end)) - mi_ini * 60 - round(sec_ini)) / 60;
            }

            if (l1cinput->verbose && gtime > gransize) 
            {
                 cout << "WARNING: rounding errors --gtime = " << gtime << " is greater than granule size = " << gransize << endl;                         }
            if ((gtime % gransize) == 0)
                gfull = "1";
            else
                gfull = "0";
             if(l1cinput->verbose)cout << "gtime.." << gtime << "gfull.." << gfull << endl;

            // # gridlines per granule
            // gridline index 0-num_gridlines
            ngridlines = gdindex[gran][1] - gdindex[gran][0] + 1;
            totlines += ngridlines;
            if(l1cinput->verbose)cout << "ngridlines..." << ngridlines << "tot gridlines..." << totlines << endl;
            
           if(twodays)//elapsed from previous day
           {    
           h_zero=0;
           mi_zero=0;
           sec_zero=0;
           unix2ymdhms(time_zero, &y_zero, &mo_zero, &d_zero, &h_zero, &mi_zero, &sec_zero);           
           }
           else //elapsed from the same day
           {
           d_zero=d_ini;
           mo_zero=mo_ini;
           y_zero=y_ini;
           }

            daystr = std::to_string(d_zero);
            if(daystr.size() == 1)
                daystr = "0" + daystr;
            monstr = std::to_string(mo_zero);
            if(monstr.size() == 1)
                monstr = "0" + monstr;
            yearstr = std::to_string(y_zero);

            datezerotime="seconds since " + yearstr + "-" + monstr + "-" + daystr;            
            string datemid=  unix2isodate(tg_mid, 'G');
            isodate2ydmsec((char *) datemid.substr(0,19).c_str(), &iyear, &iday, &msec);
            double dist_es=esdist_(&iyear,&iday,&msec);

            if(l1cinput->verbose)cout << "sun_earth_distance -- mid gridline" <<dist_es<<endl;
          
            // output file with list of granules
            string dateini=  unix2isodate(tg_ini, 'G');
            string datend=  unix2isodate(tg_end, 'G');

          
            yearstr=dateini.substr(0,4);
            monstr=dateini.substr(5,2);
            daystr=dateini.substr(8,2);
            hstr=dateini.substr(11,2);
            mistr=dateini.substr(14,2);
            secstr=dateini.substr(17,2);

            fdatetimestr1 = yearstr + monstr + daystr + "T" + hstr + mistr + secstr; 
            fname_out = pathstr + "PACE." + fdatetimestr1 + ".L1C" + extstr;
            fname_out_nopath = "PACE." + fdatetimestr1 + ".L1C" + extstr;

            if(l1cinput->verbose)cout << "granule filename.." << fname_out << endl;

            outf << fname_out_nopath << "," << dateini.substr(0,19) << "," << datend.substr(0,19) << "," << gfull << "\n";


            l1cfile->gridname = fname_out.c_str();
            filename_lt = fname_out.c_str();
            char* gridchar = strdup(filename_lt);
            string l1c_str = filename_lt;

            NcFile* nc_output;
            try {
                nc_output = new NcFile(filename_lt, NcFile::replace);
            } catch (NcException& e) {
                e.what();
                cerr << "l1cgen l1c_pflag= 5 : producing L1C grid: " + l1c_str << endl;
                exit(1);
            }

            bin_str binl1c;
            binl1c.sensor = l1cinput->sensor;
            meta_l1c_grid(gridchar, &binl1c, ngridlines, nc_output);

            if(l1cinput->pversion[0])
                nc_output->putAtt("processing_version", l1cinput->pversion);
            if(l1cinput->doi[0]) {
                nc_output->putAtt("identifier_product_doi", l1cinput->doi);
                nc_output->putAtt("identifier_product_doi_authority", "http://dx.doi.org");
            }
            nc_output->putAtt("history", l1cinput->history);
            nc_output->putAtt("product_name", fname_out_nopath);
            nc_output->putAtt("time_coverage_start", dateini.substr(0,19) + "Z");
            nc_output->putAtt("time_coverage_end", datend.substr(0,19) + "Z");
            nc_output->putAtt("sun_earth_distance", ncFloat, dist_es);

            NY = ngridlines;
            NX = nbinx;                       
            lat_out = allocate2d_float(NY, NX);
            lon_out = allocate2d_float(NY, NX);
            alt_out = allocate2d_float(NY, NX);
            time_nad_out = (double*)calloc(NY, sizeof(double));  // time of the day in seconds

            int cc = 0,ep_shift;
             
            if(tmgv[NY1]<0) ep_shift = 3600*24;
            else  ep_shift=0;
            if(l1cinput->verbose) cout<<"twodays flag : ep_shift :"<<ep_shift<<endl;
            for (int i = NY1; i < NY2 + 1; i++) {
                for (int j = 0; j < NX; j++) {
                    lat_out[cc][j] = lat_gd[i][j];
                    lon_out[cc][j] = lon_gd[i][j];
                    alt_out[cc][j] = alt_gd[i][j];            
                    time_nad_out[cc] = tmgv[i]+ep_shift;               
                }
                cc++;
            }
            if(lat_gd[NY1+1][nadir_bin_index]>lat_gd[NY1][nadir_bin_index]) {dirstr1="Ascending";l1cfile->orb_dir=1;} else {dirstr1= "Descending";l1cfile->orb_dir=0;}
            if(lat_gd[NY2][nadir_bin_index]>lat_gd[NY2-1][nadir_bin_index]) dirstr2="Ascending";else dirstr2= "Descending";
            nc_output->putAtt("startDirection", dirstr1);
            nc_output->putAtt("endDirection", dirstr2);

            // vars
            NcGroup ba_grp = nc_output->getGroup("bin_attributes");
            NcVar v1 = ba_grp.getVar("nadir_view_time");
            v1.putVar(&time_nad_out[0]);
            v1.putAtt("units",datezerotime);
            NcGroup geo_grp = nc_output->getGroup("geolocation_data");
            v1 = geo_grp.getVar("latitude");
            v1.putVar(&lat_out[0][0]);
            v1 = geo_grp.getVar("longitude");
            v1.putVar(&lon_out[0][0]);
            v1 = geo_grp.getVar("height");
            v1.putVar(&alt_out[0][0]);

            // GRING-------
            // determine the number of GCpoint indexes
            // default 6 for the swath sides + number of coordinates every 20 degrees latitude

            Ngring = 4;
            if(Ngring>0){
                float* latarr = (float*)calloc(Ngring, sizeof(float));
                float* lonarr = (float*)calloc(Ngring, sizeof(float));
            // ascending pass
            if (l1cfile->orb_dir == 1) {
                // corners--counterclockwise and ascending pass
                        latarr[0] = lat_gd[NY2][l1cfile->nbinx - 1];              
                        latarr[1] = lat_gd[NY2][0];
                        latarr[2] = lat_gd[NY1][0];               
                        latarr[3] = lat_gd[NY1][l1cfile->nbinx - 1];

                        lonarr[0] = lon_gd[NY2][l1cfile->nbinx - 1];
                        lonarr[1] = lon_gd[NY2][0];
                        lonarr[2] = lon_gd[NY1][0];
                        lonarr[3] = lon_gd[NY1][l1cfile->nbinx - 1];
        
                } else  // descending orb
                {               
                        latarr[0] = lat_gd[NY2][l1cfile->nbinx - 1];
                        latarr[1] = lat_gd[NY2][0];
                        latarr[2] = lat_gd[NY1][0];
                        latarr[3] = lat_gd[NY1][l1cfile->nbinx - 1];

                        lonarr[0] = lon_gd[NY2][l1cfile->nbinx - 1];
                        lonarr[1] = lon_gd[NY2][0];
                        lonarr[2] = lon_gd[NY1][0];
                        lonarr[3] = lon_gd[NY1][l1cfile->nbinx - 1];

               }  // end descending

                         string onelat,onelon,onecoor,firstcoor;
                //GRing params
                string gring_polygon,gring_crs,gring_latmin,gring_latmax,gring_lonmin,gring_lonmax;
                dp=4;

                for (int s = 0; s < dp; s++) {        
                    onelat=to_string(latarr[s]);
                    onelon=to_string(lonarr[s]);                        
                   if(s==0)
                    {
                        onecoor="POLYGON(("+onelat+" "+onelon+",";
                        firstcoor=onelat+" "+onelon;
                    }
                    else onecoor=onelat+" "+onelon+",";

                    gring_polygon+=onecoor;                         
                }

                
                //last coor
                gring_polygon+=firstcoor+"))";
              
                float latpmin=999,latpmax=-999,lonpmin=999,lonpmax=-999;
               for(int row=NY1;row<NY2+1;row++)
                {
                    for(int col=0;col<nbinx;col++)
                    {

                      if(lat_gd[row][col]<latpmin && lat_gd[row][col]!= BAD_FLT) latpmin=lat_gd[row][col];
                      if(lat_gd[row][col]>latpmax && lat_gd[row][col]!= BAD_FLT) latpmax=lat_gd[row][col];
                      if(lon_gd[row][col]<lonpmin && lon_gd[row][col]!= BAD_FLT) lonpmin=lon_gd[row][col];
                      if(lon_gd[row][col]>lonpmax  && lon_gd[row][col]!= BAD_FLT) lonpmax=lon_gd[row][col];
                    }}

                gring_crs="EPSG:4326";
                nc_output->putAtt("geospatial_bounds", gring_polygon);
                nc_output->putAtt("geospatial_bounds_crs", gring_crs);

                //algo to find min/max coordinates inside the gring polygon
                nc_output->putAtt("geospatial_lat_min", ncFloat, latpmin);
                nc_output->putAtt("geospatial_lat_max", ncFloat, latpmax);
                nc_output->putAtt("geospatial_lon_min", ncFloat, lonpmin);
                nc_output->putAtt("geospatial_lon_max", ncFloat, lonpmax);
  
                free (latarr);
                free (lonarr);                                
            } else {
                cout << "ERROR EXTRACTING GRING coordinates!!-----" << endl;
                exit(1);
            }

            nc_output->close();

            if (lat_out != nullptr)
                free (lat_out);
            if (lon_out != nullptr)
                free (lon_out);
            if (alt_out != nullptr)
                free (alt_out);
            if (time_nad_out != nullptr)
                free (time_nad_out);

           add_proc_group_l1c(l1cinput,l1cfile,fname_out.c_str());  
        }  // gdtime >0
    }      // end granule loop

    if (gdtime != nullptr)
        free (gdtime);
    gdtime = nullptr;

    if (gdindex != nullptr)
        free (gdindex);
    gdindex = nullptr;

    if (granid != nullptr)
        free (granid);
    granid = nullptr;

    outf.close();

    return 0;
}




// this version includes info coming from Don open file routines--

int32_t L1C::load_l1c_filehandle4(l1c_filehandle* l1cfile, L1C_input* l1cinput) {
    int nfiles;
    string filename;
    nfiles = l1cinput->files.size();

    for (int j = 0; j < nfiles; j++) {
        filename = l1cinput->files[j];
        l1cfile->ifiles.push_back(filename);
    }

    l1cfile->l1c_pflag = l1cinput->l1c_pflag;
    l1cfile->swath_num = l1cinput->swath_num;

    for (int j = 0; j < 10; j++) {  // up to 10 granules
        l1cfile->selgran[j] = l1cinput->selgran[j];
    }
    // projection params
    l1cfile->gres = l1cinput->grid_resolution;
    l1cfile->verbose = l1cinput->verbose;
    l1cfile->proj_type = l1cinput->projection;
    l1cfile->cloud_corrected = l1cinput->cloud_correct;
    // multi  attributes (view, pol, bands)
    l1cfile->overlap_vflag = l1cinput->overlap_vflag;
    l1cfile->overlap_pflag = l1cinput->overlap_pflag;
    l1cfile->overlap_bflag = l1cinput->overlap_bflag;
    // uncertainty params l1c merged products
    l1cfile->unc_meth = l1cinput->unc_meth;
    l1cfile->unc_thres_v = l1cinput->unc_thres_v;
    l1cfile->unc_thres_p = l1cinput->unc_thres_p;
    l1cfile->unc_thres_b = l1cinput->unc_thres_b;

    // calibration attributes
    // l1cfile->Fobar=file->Fobar;//nc

    cout << "ok transfering filehandle to l1c_filehandle info.." << endl;
    return 0;
}

}
