//**********************************************************
// library for searching row/col from lat/lon of L1C file
//  Originally created by Martin Montes on 11/4/2022
//  last version 3/13/2023
//**********************************************************
#include <allocate2d.h>
#include <allocate3d.h>
#include "allocate4d.h"
#include <iostream>
#include <chrono>
#include <sys/stat.h>
#include <genutils.h>
#include <iostream>
#include "l1c_latlongrid.h"
#include <netcdf>
#include <GeographicLib/Geodesic.hpp>
#include <GeographicLib/Constants.hpp>
#include <GeographicLib/Rhumb.hpp>

#include <fstream>
#include <string>

using namespace std;
using namespace netCDF;
using namespace netCDF::exceptions;
using namespace GeographicLib;

bin_str::bin_str() {
    // global attrs
    version = "";
    history = "";
    pversion = "";
    doi = "";
    sensor = -1;

    // dimensions
    nviews = -1;
    nbands = -1;
    num_gridlines = -1;
    nbinx = -1;
    num_pixels = -1;
    nscans = -1;
    verbose = 0;
    // geo
    gdshift = -1;
    lat_gd = nullptr;
    lon_gd = nullptr;
    time_l1b = nullptr;
    time_gd = nullptr;
    alt = nullptr;

    // binned vars

    nrec_2D = nullptr;  // row/col
    nrec_3D = nullptr;
    nrec_3D_view = nullptr;
    nrec_4D_band = nullptr;  // row/col/view
    alt_mean = nullptr;
    alt_rmse = nullptr;
    alt_2D = nullptr;
    alt_diff2 = nullptr;  // row/col
    i_diff2 = nullptr;
    sca_3D = nullptr;     // row/col/view
    rot_angle = nullptr;
    QC_bitwise_4D = nullptr;  // row/col/view/bands
    QC_4D = nullptr;          // row/col/view
    I_4D = nullptr;           // row/col/view/bands
    I_noise_4D = nullptr;     // #pixels
 
    // OCIS line by line
    obs_per_view = nullptr;
    QC = nullptr;
    I = nullptr;
    I_noise = nullptr;

    l1file = nullptr;
    full_l1cgrid = "";
    l1cgrid = "";

    fillval1 = BAD_FLT;
    fillval2 = float(BAD_FLT);
    fillval3 = double(BAD_FLT);

    inpix = 0;
    outpix = 0;
    badgeo = 0;

    tini_l1c = -1;
    tend_l1c = -1;
    tini_l1b = -1;
    tend_l1b = -1;
    date_mid_grid = "";

    bintype = 0;
    outlist[0] = '\0';
    l1c_anc[0] = '\0';
    cloudem_flag = 0;
    cloud_type = 0;  // 0 water, 1 ice
    dem_flag=0;//dem 0 geoid height or L1C grid height between ellipsoid and geoid wgs84 in m, 1 orthometric height or dem
}

bin_str::~bin_str() {
}

int bin_str::close_bin(bin_str *binl1c) {
    if (binl1c->time_l1b != nullptr)
        free3d_double(binl1c->time_l1b);
    binl1c->time_l1b = nullptr;
    if (binl1c->time_gd != nullptr)
        free(binl1c->time_gd);
    binl1c->time_gd = nullptr;
    if (binl1c->alt != nullptr)
        free2d_float(binl1c->alt);
    binl1c->alt = nullptr;
    if (binl1c->lat_gd != nullptr)
        free2d_float(binl1c->lat_gd);
    binl1c->lat_gd = nullptr;
    if (binl1c->lon_gd != nullptr)
        free2d_float(binl1c->lon_gd);
    binl1c->lon_gd = nullptr;
    if (binl1c->alt_rmse != nullptr)
        free2d_float(binl1c->alt_rmse);
    binl1c->alt_rmse = nullptr;
    if (binl1c->alt_mean != nullptr)
        free2d_float(binl1c->alt_mean);
    binl1c->alt_mean = nullptr;
    if (binl1c->alt_2D != nullptr)
        free2d_short(binl1c->alt_2D);
    binl1c->alt_2D = nullptr;
    if (binl1c->alt_diff2 != nullptr)
        free2d_short(binl1c->alt_diff2);
    binl1c->alt_diff2 = nullptr;
    if (binl1c->i_diff2 != nullptr)
        free4d_float(binl1c->i_diff2);
    binl1c->i_diff2 = nullptr;
    if (binl1c->suna_3D != nullptr)
        free3d_float(binl1c->suna_3D);
    binl1c->suna_3D = nullptr;
    if (binl1c->sunz_3D != nullptr)
        free3d_float(binl1c->sunz_3D);
    binl1c->sunz_3D = nullptr;
    if (binl1c->sena_3D != nullptr)
        free3d_float(binl1c->sena_3D);
    binl1c->sena_3D = nullptr;
    if (binl1c->senz_3D != nullptr)
        free3d_float(binl1c->senz_3D);
    binl1c->senz_3D = nullptr;
    if (binl1c->sca_3D != nullptr)
        free3d_float(binl1c->sca_3D);
    binl1c->sca_3D = nullptr;
    if (binl1c->rot_angle != nullptr)
        free3d_float(binl1c->rot_angle);
    binl1c->rot_angle = nullptr;
    if (binl1c->nrec_2D != nullptr)
        free2d_short(binl1c->nrec_2D);
    binl1c->nrec_2D = nullptr;
    if (binl1c->nrec_3D != nullptr)
        free3d_short(binl1c->nrec_3D);
    binl1c->nrec_3D = nullptr;
    if (binl1c->nrec_3D_view != nullptr)
        free3d_short(binl1c->nrec_3D_view);
    binl1c->nrec_3D_view = nullptr;
    if (binl1c->nrec_4D_band != nullptr)
        free4d_float(binl1c->nrec_4D_band);
    binl1c->nrec_4D_band = nullptr;
    if (binl1c->I_4D != nullptr)
        free4d_float(binl1c->I_4D);
    binl1c->I_4D = nullptr;
    if (binl1c->I_noise_4D != nullptr)
        free4d_float(binl1c->I_noise_4D);
    binl1c->I_noise_4D = nullptr;

    return 0;
}

int bin_str::alloc_bin(bin_str *binl1c) {
    if (binl1c->verbose) {
        cout << "allocating/init high order dimensional arrays for L1C binning----------------------------"
             << endl;
        cout << "nviews.." << binl1c->nviews << "nbands.." << binl1c->nbands << "#gridlines.."
             << binl1c->num_gridlines << "nbinx.." << binl1c->nbinx << endl;
    }

    binl1c->nrec_2D = allocate2d_short(binl1c->num_gridlines, binl1c->nbinx);
    binl1c->alt = allocate2d_float(binl1c->num_gridlines, binl1c->nbinx);
    binl1c->alt_rmse = allocate2d_float(binl1c->num_gridlines, binl1c->nbinx);
    binl1c->alt_mean = allocate2d_float(binl1c->num_gridlines, binl1c->nbinx);
    binl1c->alt_2D = allocate2d_short(binl1c->num_gridlines, binl1c->nbinx);
    binl1c->alt_diff2 = allocate2d_short(binl1c->num_gridlines, binl1c->nbinx);
    if (binl1c->bintype == 0) {
        binl1c->lat_gd = allocate2d_float(binl1c->num_gridlines, binl1c->nbinx);
        binl1c->lon_gd = allocate2d_float(binl1c->num_gridlines, binl1c->nbinx);
    }

    binl1c->time_gd = (double *)calloc(binl1c->num_gridlines, sizeof(double));
    binl1c->time_l1b = allocate3d_double(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->suna_3D = allocate3d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->sunz_3D = allocate3d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->sena_3D = allocate3d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->senz_3D = allocate3d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->sca_3D = allocate3d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->rot_angle = allocate3d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->nrec_3D = allocate3d_short(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);
    binl1c->nrec_3D_view = allocate3d_short(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);

    binl1c->nrec_4D_band =
        allocate4d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews, binl1c->nbands);
    binl1c->I_4D = allocate4d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews, binl1c->nbands);
    binl1c->I_noise_4D =
        allocate4d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews, binl1c->nbands);
    binl1c->i_diff2 = allocate4d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews, binl1c->nbands);

    // init
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        binl1c->time_gd[i] = 0;
        for (int j = 0; j < binl1c->nbinx; j++) {
            binl1c->nrec_2D[i][j] = 0;
            binl1c->alt[i][j] = 0.;
            binl1c->alt_rmse[i][j] = 0.;
            binl1c->alt_mean[i][j] = 0.;
            binl1c->alt_2D[i][j] = 0;
            binl1c->alt_diff2[i][j] = 0;
            if (binl1c->bintype == 0) {
                binl1c->lat_gd[i][j] = 0.;
                binl1c->lon_gd[i][j] = 0.;
            }

            for (int v = 0; v < binl1c->nviews; v++) {
                binl1c->time_l1b[i][j][v] = 0.;
                binl1c->suna_3D[i][j][v] = 0.;
                binl1c->sunz_3D[i][j][v] = 0.;
                binl1c->sena_3D[i][j][v] = 0.;
                binl1c->senz_3D[i][j][v] = 0.;
                binl1c->sca_3D[i][j][v] = 0.;
                binl1c->rot_angle[i][j][v] = 0.;
                binl1c->nrec_3D[i][j][v] = 0;
                binl1c->nrec_3D_view[i][j][v] = 0;
                for (int sb = 0; sb < binl1c->nbands; sb++) {
                    binl1c->nrec_4D_band[i][j][v][sb] = 0.;
                    binl1c->I_4D[i][j][v][sb] = 0.;
                    binl1c->I_noise_4D[i][j][v][sb] = 0.;
                    binl1c->i_diff2[i][j][v][sb] = 0;
                }
            }
        }
    }

    return 0;
}

int rmse_l1c_alt(filehandle *l1file, bin_str *binl1c, l1str *l1rec, short **gdindex) {
    short gd_row = -1, gd_col = -1;
    float fill_tilt = BAD_FLT;
    short fill_height = BAD_FLT;
    int npix = l1file->npix;

    for (int pix = 0; pix < npix; pix++) 
   {
        gd_row = gdindex[pix][0] - 1;
        gd_col = gdindex[pix][1] - 1;

        if (l1rec->height[pix] != fill_height && l1rec->tilt != fill_tilt && gd_row >= 0 && gd_col >= 0 &&
            l1rec->lat[pix] != binl1c->fillval2 && l1rec->lon[pix] != binl1c->fillval2)
       {

            binl1c->alt_diff2[gd_row][gd_col] =
                binl1c->alt_diff2[gd_row][gd_col] + (l1rec->height[pix] - binl1c->alt[gd_row][gd_col]) *
                                                        (l1rec->height[pix] - binl1c->alt[gd_row][gd_col]);        
       }   
   }

    return 0;
}

int meta_l1c_altvar(bin_str *binl1c, NcFile *nc_output) {
    float term;

    if (binl1c->verbose)
        cout << "adding alt_rmse to L1C file..." << endl;

    NcGroup geo_grp = nc_output->getGroup("geolocation_data");
    NcVar v1 = geo_grp.getVar("height_stdev");

    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            term = binl1c->alt_diff2[i][j];
            if (term >= 0 && binl1c->nrec_2D[i][j] > 0) {
                binl1c->alt_rmse[i][j] = sqrt(binl1c->alt_diff2[i][j] / binl1c->nrec_2D[i][j]);
            } else
                binl1c->alt_rmse[i][j] = binl1c->fillval2;
        }
    }
    v1.putVar(&binl1c->alt_rmse[0][0]);


 //i_stdev
   float ****temp2 = allocate4d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews, binl1c->nbands);

    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                for (int sb = 0; sb < binl1c->nbands; sb++) {
                    if (binl1c->nrec_4D_band[i][j][v][sb] > 0) { 
                        binl1c->i_diff2[i][j][v][sb]-=0.01;                                     
                        temp2[i][j][v][sb] = sqrt(binl1c->i_diff2[i][j][v][sb]/binl1c->nrec_4D_band[i][j][v][sb]);
                    } else {
                        temp2[i][j][v][sb] = binl1c->fillval2;
                    }
                }
            }
        }
    }

    NcGroup od_grp = nc_output->getGroup("observation_data");
    v1 = od_grp.getVar("i_stdev"); 
    if (binl1c->bintype == 0) {
        v1 = od_grp.getVar("i_stdev");
        v1.putVar(&temp2[0][0][0][0]);
    }

    free4d_float(temp2);

    return 0;
}

int meta_l1c_bin(filehandle *l1file, bin_str *binl1c, NcFile *nc_output) {
    if (binl1c->verbose)
        cout << "adding final binned vars to ...." << binl1c->full_l1cgrid << endl;

    NcDim ydim = nc_output->getDim("bins_along_track");
    NcDim xdim = nc_output->getDim("bins_across_track");
    NcDim vdim = nc_output->getDim("number_of_views");

    std::vector<NcDim> dimvec3;
    dimvec3.push_back(ydim);
    dimvec3.push_back(xdim);
    dimvec3.push_back(vdim);

    NcGroup geo_grp = nc_output->getGroup("geolocation_data");

    // time_offsets

    double ***time_offsets = nullptr, toff, toff_fill, toff_min, toff_max;


    time_offsets = allocate3d_double(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);

    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                time_offsets[i][j][v] = 0;
            }
        }
    }

    NcGroup ba_grp = nc_output->getGroup("bin_attributes");
    NcVar v1 = ba_grp.getVar("nadir_view_time");
    v1.getVar(&binl1c->time_gd[0]);

    v1 = ba_grp.getVar("view_time_offsets");
    NcVarAtt a1 = v1.getAtt("_FillValue");  // root group
    a1.getValues(&toff_fill);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&toff_min);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&toff_max);

    // sensor azimuth---
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                if (binl1c->nrec_3D[i][j][v] > 0) {
                    toff = binl1c->time_gd[i] - (binl1c->time_l1b[i][j][v] / binl1c->nrec_3D[i][j][v]);

                    if (toff >= toff_min && toff <= toff_max) {
                        time_offsets[i][j][v] = toff;
                    } else
                        time_offsets[i][j][v] = toff_fill;
                } else {
                    time_offsets[i][j][v] = toff_fill;
                }
            }
        }
    }

    v1 = ba_grp.getVar("view_time_offsets");
    v1.putVar(&time_offsets[0][0][0]);

    free3d_double(time_offsets);

    v1 = geo_grp.getVar("sensor_azimuth_angle");

    short ***temp;
    float angleScale = 100.0;

    temp = allocate3d_short(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews);

    // sensor azimuth---
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                if (binl1c->nrec_3D[i][j][v] > 0) {
                    temp[i][j][v] = (binl1c->sena_3D[i][j][v] / binl1c->nrec_3D[i][j][v]) * angleScale;
                } else {
                    temp[i][j][v] = binl1c->fillval1;
                }
            }
        }
    }

    v1.putVar(&temp[0][0][0]);

    // sensor zenith---
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                if (binl1c->nrec_3D[i][j][v] > 0) {
                    temp[i][j][v] = (binl1c->senz_3D[i][j][v] / binl1c->nrec_3D[i][j][v] )* angleScale;
                } else {
                    temp[i][j][v] = binl1c->fillval1;
                }
            }
        }
    }

    v1 = geo_grp.getVar("sensor_zenith_angle");
    v1.putVar(&temp[0][0][0]);

    // solar azimuth---
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                if (binl1c->nrec_3D[i][j][v] > 0) {
                    temp[i][j][v] = (binl1c->suna_3D[i][j][v] / binl1c->nrec_3D[i][j][v] )* angleScale;
                } else {
                    temp[i][j][v] = binl1c->fillval1;
                }
            }
        }
    }

    v1 = geo_grp.getVar("solar_azimuth_angle");
    v1.putVar(&temp[0][0][0]);

    // sensor zenith---
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                if (binl1c->nrec_3D[i][j][v] > 0) {
                    temp[i][j][v] = (binl1c->sunz_3D[i][j][v] / binl1c->nrec_3D[i][j][v] )* angleScale;
                } else {
                    temp[i][j][v] = binl1c->fillval1;
                }
            }
        }
    }

    v1 = geo_grp.getVar("solar_zenith_angle");
    v1.putVar(&temp[0][0][0]);

    // scattering angle---
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                if (binl1c->nrec_3D[i][j][v] > 0) {
                    temp[i][j][v] = (binl1c->sca_3D[i][j][v] / binl1c->nrec_3D[i][j][v]) * angleScale;
                } else {
                    temp[i][j][v] = binl1c->fillval1;
                }
            }
        }
    }

    v1 = geo_grp.getVar("scattering_angle");
    v1.putVar(&temp[0][0][0]);

    // rotation angle---
    // for (int i = 0; i < binl1c->num_gridlines; i++) {
    //     for (int j = 0; j < binl1c->nbinx; j++) {
    //         for (int v = 0; v < binl1c->nviews; v++) {
    //             if (binl1c->nrec_3D[i][j][v] > 0) {
    //                 temp[i][j][v] = (binl1c->rot_angle[i][j][v] / binl1c->nrec_3D[i][j][v]) * angleScale;
    //             } else {
    //                 temp[i][j][v] = binl1c->fillval1;
    //             }
    //         }
    //     }
    // }

    // v1 = geo_grp.getVar("rotation_angle");
    // v1.putVar(&temp[0][0][0]);
    free3d_short(temp);

    NcGroup od_grp = nc_output->getGroup("observation_data");
    v1 = od_grp.getVar("number_of_observations");  // ONLY 1 BAND!!!!!! NO QUALITY CONTROL,constrained no fillvalues,
                                         // and withing min/max Lt

    if (binl1c->bintype == 0) {
        v1.putVar(&binl1c->nrec_3D[0][0][0]);
    }

    // I/reflectance
    float ****temp2 = allocate4d_float(binl1c->num_gridlines, binl1c->nbinx, binl1c->nviews, binl1c->nbands);
    for (int i = 0; i < binl1c->num_gridlines; i++) {
        for (int j = 0; j < binl1c->nbinx; j++) {
            for (int v = 0; v < binl1c->nviews; v++) {
                for (int sb = 0; sb < binl1c->nbands; sb++) {
                    if (binl1c->nrec_4D_band[i][j][v][sb] > 0) {
                        temp2[i][j][v][sb] = binl1c->I_4D[i][j][v][sb] / binl1c->nrec_4D_band[i][j][v][sb];                        
//                                                cout<<"v "<<v+1<<"sb "<<sb+1<<"row "<<i+1<<"col.."<<j+1<<"MEAN radiance.. OCI = "<<temp2[i][j][v][sb]<<endl;                                                                   
                    } else {
                        temp2[i][j][v][sb] = binl1c->fillval2;
                    }
                }
            }
        }
    }

    if (binl1c->bintype == 0) {
        v1 = od_grp.getVar("i");
        v1.putVar(&temp2[0][0][0][0]);
    }

    free4d_float(temp2);

    return 0;
}

int meta_l1c_full(filehandle *l1file, bin_str *binl1c, const char *l1c_grid, NcFile *nc_output) {
    string senstr, titlestr, prodstr, ifile_str;
    int NVIEWS, NBANDS;  // NBANDS_POL;
    int32_t xbins;       //,ybins;
    char *ifile_char = l1file->name;
    file_format format = getFormat(ifile_char);
    int16_t syear, smon, sday, syear2, smon2, sday2;
    double secs, secs2;
    double tgran_ini, tgran_ini_sec, tgran_end, tgran_end_sec;
    string gfull;
    int16_t gtime;
    int16_t y_ini, mo_ini, d_ini, h_ini, mi_ini, y_end, mo_end, d_end, h_end, mi_end;
    double sec_ini, sec_end;
    int32_t gransize = 5;
    int logoff = -1;

    if (binl1c->outlist[0] == '\0')
        strcpy(binl1c->outlist, "l1c.tmp");
    string outxt(binl1c->outlist);

    std::ofstream outf;

    outf.open(outxt, std::ofstream::out | std::ofstream::trunc);

    if (outf) {
        if (binl1c->verbose)
            cout << "writing L1C granules to outfile..." << outxt << endl;
    } else {
        std::cerr << "output file.." << outxt << " could not be opened for writing!\n";
        return 1;
    }

    if (format.type == FT_HKT || format.type == FT_L1C)
        format.type = FT_OCIL1B;  // if HKT then FT_OCIS----

    if (format.type == FT_SPEXONE) {
        senstr = "SPEXONE";
        titlestr = "PACE SPEXone Level-1C Data";
        xbins = 29;
        NVIEWS = 5;
        NBANDS = 400;
        //    NBANDS_POL=50;
    } else if (format.type == FT_HARP2) {
        senstr = "HARP2";
        titlestr = "PACE HARP2 Level-1C Data";
        xbins = 457;
        NVIEWS = 90;
        NBANDS = 1;
        //    NBANDS_POL=1;
    } else if (format.type == FT_OCIL1B) {
        senstr = "OCI";
        titlestr = "PACE OCI Level-1C Data";
        xbins = 519;
        NVIEWS = 2;
        NBANDS = 286;
    } else  // OCI
    {
        senstr = "OCI";
        titlestr = "PACE OCI Level-1C Data";
        xbins = 519;
        NVIEWS = 2;
        NBANDS = 286;
    }

    // l1c grid
    string l1c_str = l1c_grid;

    NcFile *nc_l1cgrid;
    try {
        nc_l1cgrid = new NcFile(l1c_grid, NcFile::read);
    } catch (NcException &e) {
        e.what();
        cerr << "l1cgen l1c_pflag= 8:: Failure reading L1C grid: " + l1c_str << endl;
        exit(1);
    }

    // l1c_grid--------
    NcDim yd = nc_l1cgrid->getDim("bins_along_track");
    int32_t num_gridlines = yd.getSize();
    binl1c->nviews = NVIEWS;
    binl1c->nbands = NBANDS;
    binl1c->num_gridlines = num_gridlines;
    binl1c->nbinx = xbins;

    // alloc binned vars in binl1c class
    binl1c->alloc_bin(binl1c);

    NcDim vdim = nc_l1cgrid->getDim("number_of_views");

    //   NcDim idim=nc_l1cgrid->getDim("intensity_bands_per_view");

    meta_l1c_global(ifile_char, binl1c,num_gridlines, nc_output);

    NcDim idim = nc_output->getDim("intensity_bands_per_view");

    // adding fill values of binned variables----
    // sena,senz,sola,solz,sca, I
    NcGroup geo_grp = nc_output->getGroup("geolocation_data");
    NcVar v1 = geo_grp.getVar("sensor_azimuth_angle");
    short value;
    NcVarAtt a1 = v1.getAtt("_FillValue");  // root group
    a1.getValues(&value);

    binl1c->fillval1 = value;  //-32767s

    NcGroup od_grp = nc_output->getGroup("observation_data");
    v1 = od_grp.getVar("i");
    float value2;
    a1 = v1.getAtt("_FillValue");  // root group
    a1.getValues(&value2);

    binl1c->fillval2 = value2;

    // global attributes
    if(!binl1c->pversion.empty())
        nc_output->putAtt("processing_version", binl1c->pversion);
    if(!binl1c->doi.empty()) {
        nc_output->putAtt("identifier_product_doi", binl1c->doi);
        nc_output->putAtt("identifier_product_doi_authority", "http://dx.doi.org");
    }

    nc_output->putAtt("history", binl1c->history);
    string name;
    nc_output->putAtt("product_name", binl1c->full_l1cgrid);
    NcGroupAtt i1 = nc_l1cgrid->getAtt("startDirection");  // NcGroupAtt is a global attr!!
    i1.getValues(name);
    nc_output->putAtt("startDirection", name);
    i1 = nc_l1cgrid->getAtt("endDirection");
    i1.getValues(name);
    nc_output->putAtt("endDirection", name);
    i1 = nc_l1cgrid->getAtt("time_coverage_start");
    i1.getValues(name);
    nc_output->putAtt("time_coverage_start", name);
    tgran_ini = isodate2unix(name.c_str());
    i1 = nc_l1cgrid->getAtt("time_coverage_end");
    i1.getValues(name);
    tgran_end = isodate2unix(name.c_str());
    nc_output->putAtt("time_coverage_end", name);
    i1 = nc_l1cgrid->getAtt("sun_earth_distance");
    i1.getValues(&value2);
    nc_output->putAtt("sun_earth_distance", ncFloat, value2);

    unix2ymds(tgran_ini, &syear, &smon, &sday, &secs);
    unix2ymds(tgran_end, &syear2, &smon2, &sday2, &secs2);
    tgran_ini_sec = ymds2unix(syear, smon, sday, secs);
    tgran_end_sec = ymds2unix(syear2, smon2, sday2, secs2);

    unix2ymdhms(tgran_ini_sec, &y_ini, &mo_ini, &d_ini, &h_ini, &mi_ini, &sec_ini);
    unix2ymdhms(tgran_end_sec, &y_end, &mo_end, &d_end, &h_end, &mi_end, &sec_end);

    if (mi_end * 60 == 0) {
        gtime = ((60 * 60 + round(sec_end)) - mi_ini * 60 - round(sec_ini)) / 60;
    } else {
        gtime = ((mi_end * 60 + round(sec_end)) - mi_ini * 60 - round(sec_ini)) / 60;
    }

    if (gtime > gransize) {
        cout << "gtime = " << gtime << " is greater than granule size = " << gransize << endl;
        exit(1);
    }
    if ((gtime % gransize) == 0)
        gfull = "1";
    else
        gfull = "0";

    // time coverage start----
    string secstr = std::to_string(sec_ini);
    string mistr = std::to_string(mi_ini);
    string hstr = std::to_string(h_ini);
    string daystr = std::to_string(d_ini);
    string monstr = std::to_string(mo_ini);
    string yearstr = std::to_string(y_ini);

    int length = (int)floor(log10(mo_ini)) + 1;
    if (length == 1)
        monstr = "0" + monstr;
    length = (int)floor(log10(d_ini)) + 1;
    if (length == 1)
        daystr = "0" + daystr;

    if (h_ini == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(h_ini + logoff)) + 1;
    if (length == 1)
        hstr = "0" + hstr;
    if (mi_ini == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(mi_ini + logoff)) + 1;
    if (length == 1)
        mistr = "0" + mistr;
    if (sec_ini == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(round(sec_ini + logoff))) + 1;
    if (length == 1)
        secstr = "0" + secstr;

    string fdatetimestr1 = yearstr + monstr + daystr + "T" + hstr + mistr + secstr.substr(0, 2);
    string datetimestr1 =
        yearstr + "-" + monstr + "-" + daystr + "T" + hstr + ":" + mistr + ":" + secstr.substr(0, 2);
    // time coevrage end----
 
    secstr = std::to_string(sec_end);
    mistr = std::to_string(mi_end);
    hstr = std::to_string(h_end);
    daystr = std::to_string(d_end);
    monstr = std::to_string(mo_end);
    yearstr = std::to_string(y_end);

    length = (int)floor(log10(mo_end)) + 1;
    if (length == 1)
        monstr = "0" + monstr;
    length = (int)floor(log10(d_end)) + 1;
    if (length == 1)
        daystr = "0" + daystr;

    if (h_end == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(h_end + logoff)) + 1;
    if (length == 1)
        hstr = "0" + hstr;
    if (mi_end == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(mi_end + logoff)) + 1;
    if (length == 1)
        mistr = "0" + mistr;
    if (sec_end == 0)
        logoff = 1;
    else
        logoff = 0;
    length = (int)floor(log10(round(sec_end + logoff))) + 1;
    if (length == 1)
        secstr = "0" + secstr;

    string datetimestr2 =
        yearstr + "-" + monstr + "-" + daystr + "T" + hstr + ":" + mistr + ":" + secstr.substr(0, 2);
    string extstr = ".nc";
    string fname_out_nopath = binl1c->full_l1cgrid;
    outf << fname_out_nopath << "," << datetimestr1 << "," << datetimestr2 << "," << gfull << "\n";
    outf.close();

    double *time_nad = (double *)calloc(num_gridlines, sizeof(double));
    NcGroup ba_grp = nc_l1cgrid->getGroup("bin_attributes");
    v1 = ba_grp.getVar("nadir_view_time");
    v1.getVar(time_nad);
    string tmpUnits;
    v1.getAtt("units").getValues(tmpUnits);
    ba_grp = nc_output->getGroup("bin_attributes");
    v1 = ba_grp.getVar("nadir_view_time");
    v1.putVar(&time_nad[0]);
    if(!tmpUnits.empty())
        v1.putAtt("units", tmpUnits);
    if (time_nad != nullptr)
        free(time_nad);
    time_nad = nullptr;

    geo_grp = nc_l1cgrid->getGroup("geolocation_data");

    if (binl1c->bintype == 0) {
        v1 = geo_grp.getVar("latitude");
        v1.getVar(&binl1c->lat_gd[0][0]);
        geo_grp = nc_output->getGroup("geolocation_data");
        v1 = geo_grp.getVar("latitude");
        v1.putVar(&binl1c->lat_gd[0][0]);

        geo_grp = nc_l1cgrid->getGroup("geolocation_data");
        v1 = geo_grp.getVar("longitude");
        v1.getVar(&binl1c->lon_gd[0][0]);
        geo_grp = nc_output->getGroup("geolocation_data");
        v1 = geo_grp.getVar("longitude");
        v1.putVar(&binl1c->lon_gd[0][0]);
    }

    geo_grp = nc_l1cgrid->getGroup("geolocation_data");
    v1 = geo_grp.getVar("height");
    //   v1=geo_grp.getVar("altitude");
    v1.getVar(&binl1c->alt[0][0]);
    geo_grp = nc_output->getGroup("geolocation_data");
    v1 = geo_grp.getVar("height");
    v1.putVar(&binl1c->alt[0][0]);

    // GRING----
    string gs_bounds, crs;
    float gs_lat_max=-999, gs_lat_min=-999, gs_lon_max=-999, gs_lon_min=-999;
    i1 = nc_l1cgrid->getAtt("geospatial_bounds");
    i1.getValues(gs_bounds);
    nc_output->putAtt("geospatial_bounds", gs_bounds);

    i1 = nc_l1cgrid->getAtt("geospatial_bounds_crs");
    i1.getValues(crs);
    nc_output->putAtt("geospatial_bounds_crs", crs);

    i1 = nc_l1cgrid->getAtt("geospatial_lat_max");
    i1.getValues(&gs_lat_max);
    nc_output->putAtt("geospatial_lat_max", ncFloat, gs_lat_max);

    i1 = nc_l1cgrid->getAtt("geospatial_lat_min");
    i1.getValues(&gs_lat_min);
    nc_output->putAtt("geospatial_lat_min", ncFloat, gs_lat_min);

    i1 = nc_l1cgrid->getAtt("geospatial_lon_max");
    i1.getValues(&gs_lon_max);
    nc_output->putAtt("geospatial_lon_max", ncFloat, gs_lon_max);

    i1 = nc_l1cgrid->getAtt("geospatial_lon_min");
    i1.getValues(&gs_lon_min);
    nc_output->putAtt("geospatial_lon_min", ncFloat, gs_lon_min);

    std::vector<NcDim> dimvec2_rad;
    dimvec2_rad.push_back(vdim);
    dimvec2_rad.push_back(idim);

    float *Fobar = l1file->Fobar;
    float *fwave = l1file->fwave;

    float fwave_view[NVIEWS][NBANDS], fwhm_view[NVIEWS][NBANDS], Fobar_view[NVIEWS][NBANDS];

    for (int i = 0; i < NVIEWS; i++) {
        for (int j = 0; j < NBANDS; j++) {
            fwave_view[i][j] = fwave[j];
            fwhm_view[i][j] = 5.0;
            Fobar_view[i][j] = 10 * Fobar[j];
        }
    }

    // intensity_wavelengths, bandpasses and F0---
    NcGroup svb_grp = nc_output->getGroup("sensor_views_bands");
    v1 = svb_grp.addVar("intensity_wavelength", ncFloat, dimvec2_rad);  // this is infull meta
    string longName = "Intensity field center wavelengths at each view";
    v1.putAtt("long_name", longName);
    string units = "nm";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat, binl1c->fillval2);
    float valid_min = 320;
    v1.putAtt("valid_min", ncFloat, valid_min);
    float valid_max = 2260;
    v1.putAtt("valid_max", ncFloat, valid_max);
    v1.putVar(fwave_view);

    // intensity_wavelengths, bandpasses and F0---
    v1 = svb_grp.addVar("intensity_bandpass", ncFloat, dimvec2_rad);
    longName = "Intensity field bandpass at each view";
    v1.putAtt("long_name", longName);
    units = "nm";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat, binl1c->fillval2);
    valid_min = 2.5;
    v1.putAtt("valid_min", ncFloat, valid_min);
    valid_max = 100;
    v1.putAtt("valid_max", ncFloat, valid_max);
    v1.putVar(fwhm_view);

    v1 = svb_grp.addVar("intensity_f0", ncFloat, dimvec2_rad);
    longName = "Intensity band solar irradiance";
    v1.putAtt("long_name", longName);
    units = "W m^-2 um^-1";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat,binl1c->fillval2);
    valid_min = 0;
    v1.putAtt("valid_min", ncFloat, valid_min);
    valid_max = 4000;
    v1.putAtt("valid_max", ncFloat, valid_max);
    v1.putVar(Fobar_view);

    nc_l1cgrid->close();

    return 0;
}

int meta_l1c_grid(char *gridname, bin_str *binl1c, int16_t num_gridlines, NcFile *nc_output) {
    string senstr, titlestr, prodstr, ifile_str;
    string date_created;
    int NVIEWS, NBANDS, NBANDS_POL;
    int32_t xbins, ybins;
    int nadir_bin_index;
    ybins = num_gridlines;

    if (binl1c->sensor == 34) {
        senstr = "SPEXONE";
        titlestr = "PACE SPEXone Level-1C Data";
        nadir_bin_index = 14;
        xbins = 29;
        NVIEWS = 5;
        NBANDS = 400;
        NBANDS_POL = 50;
    } else if (binl1c->sensor == 35) {
        senstr = "HARP2";
        titlestr = "PACE HARP2 Level-1C Data";
        nadir_bin_index = 228;
        xbins = 457;
        NVIEWS = 90;
        NBANDS = 1;
        NBANDS_POL = 1;
    } else if (binl1c->sensor == 30 || binl1c->sensor == 31) {
        senstr = "OCI";
        titlestr = "PACE OCI Level-1C Data";
        nadir_bin_index = 259;
        xbins = 519;
        NVIEWS = 2;
        NBANDS = 239;  // 249 originally no polarization bands
    } else
    {
        senstr = "OCI";
        titlestr = "PACE OCI Level-1C Data";
        nadir_bin_index =259;
        xbins = 519;
        NVIEWS = 2;
        NBANDS = 239;  // 249 ORIGINALLY no polarization bands
    }

    prodstr = string(gridname);
    // get rid of sensor--
    if (binl1c->verbose)
        cout << "sensor.." << senstr << "xbins.." << xbins << endl;

    // creation date---
    date_created = unix2isodate(now(), 'G');

    // global attributes---

    nc_output->putAtt("title", titlestr);
    nc_output->putAtt("instrument", senstr);
    nc_output->putAtt("Conventions", "CF-1.8 ACDD-1.3");
    nc_output->putAtt("institution", "NASA Goddard Space Flight Center, Ocean Biology Processing Group");
    nc_output->putAtt("license",
                      "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/");
    nc_output->putAtt("naming_authority", "gov.nasa.gsfc.sci.oceancolor");
    nc_output->putAtt("keywords_vocabulary", "NASA Global Change Master Directory (GCMD) Science Keywords");
    nc_output->putAtt("stdname_vocabulary", "NetCDF Climate and Forecast (CF) Metadata Convention");
    nc_output->putAtt("creator_name", "NASA/GSFC");
    nc_output->putAtt("creator_email", "data@oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("creator_url", "https://oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("project", "PACE Project");
    nc_output->putAtt("publisher_name", "NASA/GSFC");
    nc_output->putAtt("publisher_email", "data@oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("publisher_url", "https://oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("processing_level", "L1C");
    nc_output->putAtt("cdm_data_type", "swath");
    nc_output->putAtt("product_name", prodstr);
    nc_output->putAtt("date_created", date_created);
    nc_output->putAtt("sun_earth_distance", ncFloat, BAD_FLT);
    nc_output->putAtt("nadir_bin", NC_INT,nadir_bin_index);
    nc_output->putAtt("bin_size_at_nadir", "5.2km2");

    NcDim ydim = nc_output->addDim("bins_along_track", ybins);
    NcDim xdim = nc_output->addDim("bins_across_track", xbins);
    NcDim vdim = nc_output->addDim("number_of_views", NVIEWS);
    NcDim idim = nc_output->addDim("intensity_bands_per_view", NBANDS);

    // dims
    std::vector<NcDim> dimvec2_geo;
    dimvec2_geo.push_back(ydim);
    dimvec2_geo.push_back(xdim);
    std::vector<NcDim> dimvec2_rad;
    dimvec2_rad.push_back(vdim);
    dimvec2_rad.push_back(idim);
    std::vector<NcDim> dimvec3;
    dimvec3.push_back(ydim);
    dimvec3.push_back(xdim);
    dimvec3.push_back(vdim);
    std::vector<NcDim> dimvec4;
    dimvec4.push_back(ydim);
    dimvec4.push_back(xdim);
    dimvec4.push_back(vdim);
    dimvec4.push_back(idim);

    // groups
    NcGroup ba_grp = nc_output->addGroup("bin_attributes");
    NcGroup geo_grp = nc_output->addGroup("geolocation_data");

    if (binl1c->sensor == 34 || binl1c->sensor == 35) {  // SPEX 34, HARP2 35

        NcDim pdim = nc_output->addDim("polarization_bands_per_view", NBANDS_POL);
        std::vector<NcDim> dimvec2b_rad;
        dimvec2b_rad.push_back(vdim);
        dimvec2b_rad.push_back(pdim);
        std::vector<NcDim> dimvec4b;
        dimvec4b.push_back(ydim);
        dimvec4b.push_back(xdim);
        dimvec4b.push_back(vdim);
        dimvec4b.push_back(pdim);
    }

    NcVar v1 = ba_grp.addVar("nadir_view_time", ncDouble, ydim); 
    string longName = "Time bin was viewed at nadir view";
    v1.putAtt("long_name", longName);
    string units = "seconds since date";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncDouble, binl1c->fillval3);
    double valid_min_d = 0;
    v1.putAtt("valid_min", ncDouble, valid_min_d);
    double valid_max_d = 172800;
    v1.putAtt("valid_max", ncDouble, valid_max_d);

    v1 = geo_grp.addVar("latitude", ncFloat, dimvec2_geo);
    longName = "Latitudes of bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees_north";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat, binl1c->fillval2);
    float valid_min = -90;
    v1.putAtt("valid_min", ncFloat, valid_min);
    float valid_max = 90;
    v1.putAtt("valid_max", ncFloat, valid_max);

    v1 = geo_grp.addVar("longitude", ncFloat, dimvec2_geo);
    longName = "Longitudes of bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees_east";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat, binl1c->fillval2);
    valid_min = -180;
    v1.putAtt("valid_min", ncFloat, valid_min);
    valid_max = 180;
    v1.putAtt("valid_max", ncFloat, valid_max);

    v1 = geo_grp.addVar("height", ncShort, dimvec2_geo);
    longName = "Altitude at bin locations";
    v1.putAtt("long_name", longName);
    units = "meters";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort, binl1c->fillval1);
    // use min/max from gebco file
  //  short valid_min_s = -10952;
    short valid_min_s = -10000;
    v1.putAtt("valid_min", ncShort, valid_min_s);
//    short valid_max_s = 8627;
    short valid_max_s = 10000;
    v1.putAtt("valid_max", ncShort, valid_max_s);

    return 0;
}

int meta_l1c_global(char *gridname, bin_str *binl1c,int16_t num_gridlines, NcFile *nc_output) {
    string senstr, GATT_VAL1, prodstr,ifile_str;
    string date_created;
    int NVIEWS, NBANDS, NBANDS_POL;
    int32_t xbins, ybins;
    int nadir_bin_index;

    ybins = num_gridlines;
    prodstr = string(gridname);

    file_format format = getFormat(gridname);
    if (format.type == FT_HKT || format.type == FT_L1C)
        format.type = FT_OCIL1B;

    if (format.type == FT_SPEXONE) {
        senstr = "SPEXONE";
        GATT_VAL1 = "PACE SPEXone Level-1C Data";
        nadir_bin_index =14;
        xbins = 29;
        NVIEWS = 5;
        NBANDS = 400;
        NBANDS_POL = 50;
    } else if (format.type == FT_HARP2) {
        senstr = "HARP2";
        GATT_VAL1 = "PACE HARP2 Level-1C Data";
        nadir_bin_index =228;
        xbins = 457;
        NVIEWS = 90;
        NBANDS = 1;
        NBANDS_POL = 1;
    } else if (format.type == FT_OCIL1B) {
        senstr = "OCI";
        GATT_VAL1 = "PACE OCI Level-1C Data";
        nadir_bin_index =259;
        xbins = 519;
        NVIEWS = 2;
        NBANDS = 286;  // 249 originally no polarization bands
    } else
    {
        senstr = "OCI";
        GATT_VAL1 = "PACE OCI Level-1C Data";
        nadir_bin_index =259;
        xbins = 519;
        NVIEWS = 2;
        NBANDS = 286;  // 249 ORIGINALLY no polarization bands
    }

    // views
    float views[NVIEWS];
    if (format.type == FT_SPEXONE) {
        views[0] = -58;
        views[1] = -20;
        views[2] = 0;
        views[3] = 20;
        views[4] = 58;
    } else if (format.type == FT_OCIL1B) {
        views[0] = -20;  // OCI
        views[1] = 20;
    } else if (format.type == FT_HARP2) {
        cout << "# views TBD.....ERROR.." << endl;
        exit(1);
    } else {             // OCI default in HKT---
        views[0] = -20;  // OCI
        views[1] = 20;
    }

    // creation date---
    date_created = unix2isodate(now(), 'G');
   
    // global attributes---

    nc_output->putAtt("title", GATT_VAL1);
    nc_output->putAtt("instrument", senstr);
    nc_output->putAtt("Conventions", "CF-1.8 ACDD-1.3");
    nc_output->putAtt("institution", "NASA Goddard Space Flight Center, Ocean Biology Processing Group");
    nc_output->putAtt("license",
                      "https://science.nasa.gov/earth-science/earth-science-data/data-information-policy/");
    nc_output->putAtt("naming_authority", "gov.nasa.gsfc.sci.oceancolor");
    nc_output->putAtt("keywords_vocabulary", "NASA Global Change Master Directory (GCMD) Science Keywords");
    nc_output->putAtt("stdname_vocabulary", "NetCDF Climate and Forecast (CF) Metadata Convention");
    nc_output->putAtt("creator_name", "NASA/GSFC");
    nc_output->putAtt("creator_email", "data@oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("creator_url", "https://oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("project", "PACE Project");
    nc_output->putAtt("publisher_name", "NASA/GSFC");
    nc_output->putAtt("publisher_email", "data@oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("publisher_url", "https://oceancolor.gsfc.nasa.gov");
    nc_output->putAtt("processing_level", "L1C");
    nc_output->putAtt("cdm_data_type", "swath");
    nc_output->putAtt("product_name", prodstr);
    nc_output->putAtt("date_created", date_created);
    nc_output->putAtt("sun_earth_distance", ncFloat, BAD_FLT);
    nc_output->putAtt("nadir_bin",NC_INT,nadir_bin_index);
    nc_output->putAtt("bin_size_at_nadir", "5.2km2");

    NcDim ydim = nc_output->addDim("bins_along_track", ybins);
    NcDim xdim = nc_output->addDim("bins_across_track", xbins);
    NcDim vdim = nc_output->addDim("number_of_views", NVIEWS);
    NcDim idim = nc_output->addDim("intensity_bands_per_view", NBANDS);

    // dims
    std::vector<NcDim> dimvec2_geo;
    dimvec2_geo.push_back(ydim);
    dimvec2_geo.push_back(xdim);
    std::vector<NcDim> dimvec2_rad;
    dimvec2_rad.push_back(vdim);
    dimvec2_rad.push_back(idim);
    std::vector<NcDim> dimvec3;
    dimvec3.push_back(ydim);
    dimvec3.push_back(xdim);
    dimvec3.push_back(vdim);
    std::vector<NcDim> dimvec4;
    dimvec4.push_back(ydim);
    dimvec4.push_back(xdim);
    dimvec4.push_back(vdim);
    dimvec4.push_back(idim);

    // groups
    NcGroup svb_grp = nc_output->addGroup("sensor_views_bands");
    NcGroup ba_grp = nc_output->addGroup("bin_attributes");
    NcGroup geo_grp = nc_output->addGroup("geolocation_data");
    NcGroup od_grp = nc_output->addGroup("observation_data");


    // vars
    NcVar v1 = svb_grp.addVar("sensor_view_angle", ncFloat, vdim);
    string longName = "Along-track view angles for sensor";
    v1.putAtt("long_name", longName);
    string units = "degrees";
    v1.putAtt("units", units);    
    v1.putAtt("_FillValue", ncFloat,binl1c->fillval2);
    float valid_min = -89;
    v1.putAtt("valid_min", ncFloat, valid_min);
    float valid_max = 89;
    v1.putAtt("valid_max", ncFloat, valid_max);
    v1.putVar(views);


    if (format.type == FT_HARP2 || format.type == FT_SPEXONE) {
        NcDim pdim = nc_output->addDim("polarization_bands_per_view", NBANDS_POL);
        std::vector<NcDim> dimvec2b_rad;
        dimvec2b_rad.push_back(vdim);
        dimvec2b_rad.push_back(pdim);
        std::vector<NcDim> dimvec4b;
        dimvec4b.push_back(ydim);
        dimvec4b.push_back(xdim);
        dimvec4b.push_back(vdim);
        dimvec4b.push_back(pdim);
    }

    v1 = ba_grp.addVar("nadir_view_time", ncDouble, ydim);
    longName = "Time bin was viewed at nadir view";
    v1.putAtt("long_name", longName);
    units = "seconds since date";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncDouble, binl1c->fillval3);
    double valid_min_d = 0;
    v1.putAtt("valid_min", ncDouble, valid_min_d);
    double valid_max_d = 172800;
    v1.putAtt("valid_max", ncDouble, valid_max_d);

    v1 = ba_grp.addVar("view_time_offsets", ncDouble, dimvec3);
    longName = "Time offsets of views from nadir view";
    v1.putAtt("long_name", longName);
    units = "seconds";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncDouble, binl1c->fillval3);
    valid_min_d = -200;
    v1.putAtt("valid_min", ncDouble, valid_min_d);
    valid_max_d = 200;
    v1.putAtt("valid_max", ncDouble, valid_max_d);

    v1 = geo_grp.addVar("latitude", ncFloat, dimvec2_geo);
    longName = "Latitudes of bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees_north";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat,binl1c->fillval2);
    valid_min = -90;
    v1.putAtt("valid_min", ncFloat, valid_min);
    valid_max = 90;
    v1.putAtt("valid_max", ncFloat, valid_max);

    v1 = geo_grp.addVar("longitude", ncFloat, dimvec2_geo);
    longName = "Longitudes of bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees_east";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncFloat, binl1c->fillval2);
    valid_min = -180;
    v1.putAtt("valid_min", ncFloat, valid_min);
    valid_max = 180;
    v1.putAtt("valid_max", ncFloat, valid_max);

    v1 = geo_grp.addVar("height", ncShort, dimvec2_geo);
    longName = "Altitude at bin locations";
    v1.putAtt("long_name", longName);
    units = "meters";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort,binl1c->fillval1);
    // use min/max from gebco file
 //   short valid_min_s = -10952;
    short valid_min_s = -10000;
    v1.putAtt("valid_min", ncShort, valid_min_s);
//    short valid_max_s = 8627;
    short valid_max_s=10000;
    v1.putAtt("valid_max", ncShort, valid_max_s);

    v1 = geo_grp.addVar("height_stdev", ncShort, dimvec2_geo);
    longName = "Standard deviation of terrain altitude within bin";
    v1.putAtt("long_name", longName);
    units = "meters";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort, binl1c->fillval1);
    valid_min_s = 0;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 1000;
    v1.putAtt("valid_max", ncShort, valid_max_s);

    v1 = geo_grp.addVar("sensor_azimuth_angle", ncShort, dimvec3);
    longName = "Sensor azimuth angles at bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort, binl1c->fillval1);
    valid_min_s = -18000;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 18000;
    v1.putAtt("valid_max", ncShort, valid_max_s);
    float scale_factor = 0.01;
    v1.putAtt("scale_factor", ncFloat, scale_factor);
    float add_offset = 0;
    v1.putAtt("add_offset", ncFloat, add_offset);

    v1 = geo_grp.addVar("sensor_zenith_angle", ncShort, dimvec3);
    longName = "Sensor zenith angles at bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort,binl1c->fillval1);
    valid_min_s = 0;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 18000;
    v1.putAtt("valid_max", ncShort, valid_max_s);
    scale_factor = 0.01;
    v1.putAtt("scale_factor", ncFloat, scale_factor);
    add_offset = 0;
    v1.putAtt("add_offset", ncFloat, add_offset);

    v1 = geo_grp.addVar("solar_azimuth_angle", ncShort, dimvec3);
    longName = "Solar azimuth angle at bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort,binl1c->fillval1);
    valid_min_s = -18000;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 18000;
    v1.putAtt("valid_max", ncShort, valid_max_s);
    scale_factor = 0.01;
    v1.putAtt("scale_factor", ncFloat, scale_factor);
    add_offset = 0;
    v1.putAtt("add_offset", ncFloat, add_offset);

    v1 = geo_grp.addVar("solar_zenith_angle", ncShort, dimvec3);
    longName = "Solar zenith angle at bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort,binl1c->fillval1);
    valid_min_s = 0;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 18000;
    v1.putAtt("valid_max", ncShort, valid_max_s);
    scale_factor = 0.01;
    v1.putAtt("scale_factor", ncFloat, scale_factor);
    add_offset = 0;
    v1.putAtt("add_offset", ncFloat, add_offset);

    v1 = geo_grp.addVar("scattering_angle", ncShort, dimvec3);
    longName = "Scattering angle at bin locations";
    v1.putAtt("long_name", longName);
    units = "degrees";
    v1.putAtt("units", units);
    v1.putAtt("_FillValue", ncShort,binl1c->fillval1);
    valid_min_s = 0;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 18000;
    v1.putAtt("valid_max", ncShort, valid_max_s);
    scale_factor = 0.01;
    v1.putAtt("scale_factor", ncFloat, scale_factor);
    add_offset = 0;
    v1.putAtt("add_offset", ncFloat, add_offset);

    // v1 = geo_grp.addVar("rotation_angle", ncShort, dimvec3);
    // longName = "Rotation angle at bin locations";
    // v1.putAtt("long_name", longName);
    // units = "degrees";
    // v1.putAtt("units", units);
    // v1.putAtt("_FillValue", ncShort, binl1c->fillval1);
    // valid_min_s = 0;
    // v1.putAtt("valid_min", ncShort, valid_min_s);
    // valid_max_s = 18000;
    // v1.putAtt("valid_max", ncShort, valid_max_s);
    // scale_factor = 0.01;
    // v1.putAtt("scale_factor", ncFloat, scale_factor);
    // add_offset = 0;
    // v1.putAtt("add_offset", ncFloat, add_offset);

    v1 = od_grp.addVar("number_of_observations", ncShort, dimvec3);
    longName = "Observations contributing to bin from each view";
    v1.putAtt("long_name", longName);
    valid_min_s = 0;
    v1.putAtt("valid_min", ncShort, valid_min_s);
    valid_max_s = 999;
    v1.putAtt("valid_max", ncShort, valid_max_s);
    v1.putAtt("coordinates", "geolocation_data/longitude geolocation_data/latitude");

    if (format.type == FT_OCIL1B  || format.type == FT_HARP2) {
        uint8_t FillValue3, valid_min2, valid_max2;

        v1 = od_grp.addVar("qc", ncUbyte, dimvec4);
        longName = "quality indicator";
        v1.putAtt("long_name", longName);
        FillValue3 = 255;
        v1.putAtt("_FillValue", ncUbyte, FillValue3);
        valid_min2 = 0;
        v1.putAtt("valid_min", ncUbyte, valid_min2);
        valid_max2 = 10;
        v1.putAtt("valid_max", ncUbyte, valid_max2);
        v1.putAtt("coordinates", "geolocation_data/longitude geolocation_data/latitude");

        v1 = od_grp.addVar("i", ncFloat, dimvec4);
        longName = "I Stokes vector component";
        v1.putAtt("long_name", longName);
        units = "W m^-2 sr^-1 um^-1";
        v1.putAtt("units", units);
        v1.putAtt("_FillValue", ncFloat, binl1c->fillval2);
        valid_min = 0;
        v1.putAtt("valid_min", ncFloat, valid_min);
        valid_max = 999;
        v1.putAtt("valid_max", ncFloat, valid_max);
        v1.putAtt("coordinates", "geolocation_data/longitude geolocation_data/latitude");

        v1 = od_grp.addVar("i_stdev", ncFloat, dimvec4);
        longName = "Random stdev of i in bin";
        v1.putAtt("long_name", longName);
        units = "W m^-2 sr^-1 um^-1";
        v1.putAtt("units", units);
        v1.putAtt("_FillValue", ncFloat,binl1c->fillval2);
        valid_min = 0;
        v1.putAtt("valid_min", ncFloat, valid_min);
        valid_max = 800;
        v1.putAtt("valid_max", ncFloat, valid_max);
        v1.putAtt("coordinates", "geolocation_data/longitude geolocation_data/latitude");

    } else if (format.type == FT_SPEXONE)  // spexone
    {
    }

    return 0;
}

int bintime_l1c(filehandle *l1file, l1str *l1rec, bin_str *binstr, short **gdindex, double scantime,
                NcFile *nc_output) {
    short gd_row = -1, gd_col = -1;
    int npix = l1file->npix;
    int view;
    double timefill, mintime, maxtime;
    float fill_tilt = binstr->fillval2;
    float tilt = binstr->fillval2;

    NcGroup ba_grp = nc_output->getGroup("bin_attributes");
    NcVar v1 = ba_grp.getVar("nadir_view_time");
    NcVarAtt a1 = v1.getAtt("_FillValue");  // root group
    a1.getValues(&timefill);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&mintime);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxtime);

    // ASSUMING same time for each line!!
    for (int pix = 0; pix < npix; pix++) {
        gd_row = gdindex[pix][0] - 1;
        gd_col = gdindex[pix][1] - 1;
        tilt = l1rec->tilt;
        view = BAD_FLT;

        // 1 view at the time for OCI!!!!
        if (l1rec->tilt != fill_tilt && l1rec->tilt < -18) {  //-19.9 previous
            view = 0;
            tilt = l1rec->tilt;
        } else if (l1rec->tilt != fill_tilt && l1rec->tilt > 18) {  // 19.9 previous
            view = 1;
            tilt = l1rec->tilt;
        }

        if (gd_row >= 2000) {
            tilt = 22;
            view = 1;
            //           if(binstr->verbose) cout<<"pix # "<<pix+1<<"pos tilt beyond L1C
            //           grid"<<tilt<<"gd_row"<<gd_row<<endl;
            tilt = fill_tilt;
        }
        if (gd_row <= -1) {
            tilt = -22;
            view = 0;
            //         if(binstr->verbose)  cout<<"pix # "<<pix+1<<"neg tilt beyond L1C
            //         grid"<<tilt<<"gd_row"<<gd_row<<endl;
            tilt = fill_tilt;
        }

        //                  if(binstr->verbose)  cout<<"pix # "<<pix+1<<"scantime "<<scantime<<"min scantime
        //                  "<<mintime<<"max scantime "<<maxtime<<"gd_row "<<gd_row<<"gd_col "<<gd_col<<"tilt
        //                  "<<tilt<<endl;
        // Cumulative values---
        if (view != BAD_FLT && tilt != fill_tilt && gd_row >= 0 && gd_col >= 0 &&
            gd_row <= binstr->num_gridlines - 1 && gd_col <= binstr->nbinx - 1 && scantime != timefill &&
            scantime >= mintime && scantime <= maxtime) {  // and row<=ybin and col<xbin, I need time counter
            binstr->time_l1b[gd_row][gd_col][view] = binstr->time_l1b[gd_row][gd_col][view] + scantime;
            //                    if(binstr->verbose)  cout<<"pix # "<<pix+1<<"binning time at row
            //                    "<<gd_row+1<<"col "<<gd_col+1<<endl;
        }

    }  // end pixel

    return 0;
}

int parallax(filehandle *l1file, const char *l1c_anc, const char *l1c_grid, l1str *l1rec, bin_str *binl1c,
             short **gdindex, NcFile *nc_output,int32_t sline,int firstcall) {
    short gd_row = -1, gd_col = -1;
    string l1c_str = l1c_grid;
    string l1c_str2 = l1c_anc;
    float *lat_new = nullptr, *lon_new = nullptr;
    float dv, Re = 6378.137, Hsat = 676.5;
    float **cth = nullptr, **height = nullptr;
    float scale_factor;
    short sena_min, sena_max, senz_min, senz_max;
    int result = 0,status;
    float look_angle,*scan_angle=nullptr;//in addition to ccd bands also there is scan angle swir
    int navGrp,scangId;
    size_t start[] = { 0, 0,};
    size_t count[] = { 1, 1 };
    start[0]=sline;
    start[1]=0;
    count[0]=1;
    count[1]=l1file->npix;

    if(firstcall==0){ scan_angle = (float*)calloc(l1file->npix,sizeof(float));firstcall=1;}
    
    

    status = nc_inq_grp_ncid(l1file->sd_id, "navigation_data", &navGrp);
        check_err(status, __LINE__, __FILE__);
    status = nc_inq_varid(navGrp, "CCD_scan_angles", &scangId);
        check_err(status, __LINE__, __FILE__);
    status = nc_get_vara_float(navGrp, scangId, start, count, scan_angle);
        check_err(status, __LINE__, __FILE__);

    NcFile *nc_l1cgrid;
    try {
        nc_l1cgrid = new NcFile(l1c_grid, NcFile::read);
    } catch (NcException &e) {
        e.what();
        cerr << "l1cgen l1c_pflag= 8:: Failure reading L1C grid: " + l1c_str << endl;
        exit(1);
    }

    NcDim yd1 = nc_l1cgrid->getDim("bins_along_track");
    NcDim xd1 = nc_l1cgrid->getDim("bins_across_track");
    // get height
    NcGroup geo_grp = nc_l1cgrid->getGroup("geolocation_data");
    NcVar v1 = geo_grp.getVar("height");
    height = allocate2d_float(yd1.getSize(), xd1.getSize());
    v1.getVar(&height[0][0]);

    // anc file
   NcFile *nc_l1canc;
   if (binl1c->cloudem_flag == 1){
    try {
        nc_l1canc = new NcFile(l1c_anc, NcFile::read);
    } catch (NcException &e) {
        e.what();
        cerr << "l1cgen l1c_pflag= 8:: Failure reading L1C-ANCfile: " + l1c_str2 << endl;
        exit(1);
    }

    // 2-d as L1C grid but different nlines, get cth_ice_cloud, cth_water_cloud vars
    NcDim yd2 = nc_l1canc->getDim("bins_along_track");
    NcDim xd2 = nc_l1canc->getDim("bins_across_track");
    cth = allocate2d_float(yd2.getSize(), xd2.getSize());
    if (binl1c->cloud_type == 0) {
        v1 = nc_l1canc->getVar("cth_water_cloud");
        v1.getVar(&cth[0][0]);
    }
    if (binl1c->cloud_type == 1) {
        v1 = nc_l1canc->getVar("cth_ice_cloud");
        v1.getVar(&cth[0][0]);
    }
   }
    // compute displacement and new lat/lon for each pixel and 1 LINE!
    // compute new lat and lon

    NcGroup geo_grp2 = nc_output->getGroup("geolocation_data");
    v1 = geo_grp2.getVar("sensor_azimuth_angle");
    NcVarAtt a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_factor);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&sena_min);
    sena_min *= scale_factor;
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&sena_max);
    sena_max *= scale_factor;
    v1 = geo_grp2.getVar("sensor_zenith_angle");
    a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_factor);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&senz_min);
    senz_min *= scale_factor;
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&senz_max);
    senz_max *= scale_factor;

    lat_new = (float *)calloc(l1file->npix, sizeof(float));
    lon_new = (float *)calloc(l1file->npix, sizeof(float));

    for (int pix = 0; pix < l1file->npix; pix++) {
        gd_row = gdindex[pix][0] - 1;
        gd_col = gdindex[pix][1] - 1;
        // Displacement vector -- wang et al. 2011 ----
        if (l1rec->senz[pix] != binl1c->fillval1 && l1rec->sena[pix] != binl1c->fillval1 &&
            l1rec->senz[pix] >= senz_min && l1rec->senz[pix] <= senz_max && l1rec->sena[pix] >= sena_min &&
            l1rec->sena[pix] <= sena_max) {
            look_angle = acos(cos(l1rec->tilt * OEL_DEGRAD) * cos(scan_angle[pix] * OEL_DEGRAD));
            if (binl1c->cloudem_flag == 1)  // CTH constant in km
            {
                dv = Hsat * cth[gd_row][gd_col] * tan(look_angle) / (Hsat - cth[gd_row][gd_col]);
            } else if (binl1c->cloudem_flag == 2)  // DEM constant in km
            {
                dv =
                    Hsat * 0.001 * l1rec->height[pix] * tan(look_angle) / (Hsat - 0.001 * l1rec->height[pix]);
            }

            lat_new[pix] = l1rec->lat[pix] * OEL_DEGRAD - dv * cos(l1rec->sena[pix] * OEL_DEGRAD + OEL_PI) / Re;
            lon_new[pix] = l1rec->lon[pix] * OEL_DEGRAD - dv * sin((l1rec->sena[pix] * OEL_DEGRAD + OEL_PI)) /
                                                              (Re * cos(lat_new[pix] * OEL_DEGRAD));

            if (lon_new[pix] * OEL_RADEG > 180) {
                lon_new[pix] -= 2 * OEL_PI;
            }
            if (lon_new[pix] * OEL_RADEG < -180) {
                lon_new[pix] += 2 * OEL_PI;
            }
            lat_new[pix] *= OEL_RADEG;
            lon_new[pix] *= OEL_RADEG;
        } else {
            lat_new[pix] = binl1c->fillval2;
            lon_new[pix] = binl1c->fillval2;
        }
        if (pix == 635 && binl1c->verbose)
            cout << "sline # " << sline + 1 << "look_angle in m " << look_angle * OEL_RADEG << "dv " << dv
                 << "scan_angle " << OEL_RADEG * scan_angle[pix] << "lat_new " << lat_new[pix] << "lon_new "
                 << lon_new[pix] << "height " << l1rec->height[pix] << endl;
    }

    // recompute row/col
    result = search_l1c_parallax(l1file, lat_new, lon_new, binl1c, gdindex);

    nc_l1cgrid->close();
    if (binl1c->cloudem_flag == 1) nc_l1canc->close();
    if (cth != nullptr)
        free2d_float(cth);
    if (height != nullptr)
        free2d_float(height);
    if (lat_new != nullptr)
        free(lat_new);
    if (lon_new != nullptr)
        free(lon_new);

   if(sline==l1file->nscan-1){  if (scan_angle != nullptr) free(scan_angle);}

    return result;
}

int bin_l1c(filehandle *l1file, l1str *l1rec, bin_str *binl1c, short **gdindex, NcFile *nc_output) {
    int32_t ibp = 0, sb = 0;
    short gd_row = -1, gd_col = -1;
    int npix = l1file->npix;
    int nbands = l1file->nbands;
    int view;

    float term1, term2, term3, sca_pix, cosu, cose;
    float scale_all[7], offset_all[7];
    short minval1_all[7], maxval1_all[7];
    float minval2_all[7], maxval2_all[7];
    float fill_tilt = binl1c->fillval2;
    short fill_height = binl1c->fillval2;
    float tilt = binl1c->fillval2;
    double rotangle;
    float seed_mean;
    NcGroup geo_grp = nc_output->getGroup("geolocation_data");
    NcGroup od_grp = nc_output->getGroup("observation_data");

    // altitude
    NcVar v1 = geo_grp.getVar("height");
    scale_all[0] = 1;  // no scale and offset for height
    offset_all[0] = 0;
    NcVarAtt a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval1_all[0]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval1_all[0]);

    // sensor azimuth
    v1 = geo_grp.getVar("sensor_azimuth_angle");
    a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_all[1]);
    a1 = v1.getAtt("add_offset");  // root group
    a1.getValues(&offset_all[1]);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval1_all[1]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval1_all[1]);

    // sensor zenith
    v1 = geo_grp.getVar("sensor_zenith_angle");
    a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_all[2]);
    a1 = v1.getAtt("add_offset");  // root group
    a1.getValues(&offset_all[2]);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval1_all[2]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval1_all[2]);

    minval1_all[2] = scale_all[2] * minval1_all[2] + offset_all[2];
    maxval1_all[2] = scale_all[2] * maxval1_all[2] + offset_all[2];

    // solar zenith
    v1 = geo_grp.getVar("solar_azimuth_angle");
    a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_all[3]);
    a1 = v1.getAtt("add_offset");  // root group
    a1.getValues(&offset_all[3]);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval1_all[3]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval1_all[3]);

    // solar zenith
    v1 = geo_grp.getVar("solar_zenith_angle");
    a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_all[4]);
    a1 = v1.getAtt("add_offset");  // root group
    a1.getValues(&offset_all[4]);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval1_all[4]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval1_all[4]);

    // scattering angle
    v1 = geo_grp.getVar("scattering_angle");
    a1 = v1.getAtt("scale_factor");  // root group
    a1.getValues(&scale_all[5]);
    a1 = v1.getAtt("add_offset");  // root group
    a1.getValues(&offset_all[5]);
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval1_all[5]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval1_all[5]);

    v1 = od_grp.getVar("i");
    a1 = v1.getAtt("valid_min");  // root group
    a1.getValues(&minval2_all[0]);
    a1 = v1.getAtt("valid_max");  // root group
    a1.getValues(&maxval2_all[0]);
  
    for (int pix = 0; pix < npix; pix++) {
        gd_row = gdindex[pix][0] - 1;
        gd_col = gdindex[pix][1] - 1;
        tilt = l1rec->tilt;
        view = BAD_FLT;

        if (l1rec->tilt != fill_tilt && l1rec->tilt < -18) {  //-19.9 previous
            view = 0;
            tilt = l1rec->tilt;
        } else if (l1rec->tilt != fill_tilt && l1rec->tilt > 18) {  // 19.9 previous
            view = 1;
            tilt = l1rec->tilt;
        }
        if (gd_row >= 2000) {
            tilt = 22;
            view = 1;
            tilt = fill_tilt;
        }
        if (gd_row <= -1) {
            tilt = -22;
            view = 0;
            tilt = fill_tilt;
        }

        // 2-D arrays---needing to add fillvalue != for l1rec->alt[pix]
        if (view != BAD_FLT && tilt != fill_tilt && gd_row >= 0 && gd_col >= 0 &&
            gd_row <= binl1c->num_gridlines - 1 && gd_col <= binl1c->nbinx - 1 &&
            l1rec->lat[pix] != binl1c->fillval2 && l1rec->lon[pix] != binl1c->fillval2) {
            if (l1rec->height[pix] != fill_height)  // ADD MIN/MAX CONSTRAINT
            {
                binl1c->nrec_2D[gd_row][gd_col] += 1;
                binl1c->alt_2D[gd_row][gd_col] =
                    binl1c->alt_2D[gd_row][gd_col] + l1rec->height[pix]; 
            }
            binl1c->inpix++;
         
        } else {
            binl1c->outpix++;
        }

        // Cumulative values---
        // geometry
        if (view != BAD_FLT && tilt != fill_tilt && gd_row >= 0 && gd_col >= 0 &&
            gd_row <= binl1c->num_gridlines - 1 && gd_col <= binl1c->nbinx - 1 &&
            l1rec->senz[pix] != binl1c->fillval1 && l1rec->senz[pix] >= minval1_all[2] &&
            l1rec->senz[pix] <= maxval1_all[2]) {
            binl1c->nrec_3D[gd_row][gd_col][view] += 1;
            binl1c->sena_3D[gd_row][gd_col][view] = binl1c->sena_3D[gd_row][gd_col][view] + l1rec->sena[pix];
            binl1c->senz_3D[gd_row][gd_col][view] =
                binl1c->senz_3D[gd_row][gd_col][view] + l1rec->senz[pix];  // sensor zenith goes from 0 to 180
            binl1c->suna_3D[gd_row][gd_col][view] = binl1c->suna_3D[gd_row][gd_col][view] + l1rec->sola[pix];
            binl1c->sunz_3D[gd_row][gd_col][view] = binl1c->sunz_3D[gd_row][gd_col][view] + l1rec->solz[pix];

            cose = cos((l1rec->senz[pix] + 180) * OEL_DEGRAD);
            cosu = cos((l1rec->solz[pix]) * OEL_DEGRAD);
            term1 = cose * cosu;
            term2 = sqrt((1 - cose * cose)) * sqrt((1 - cosu * cosu));
            term3 = cos((l1rec->sena[pix] + 180 - l1rec->sola[pix]) * OEL_DEGRAD);
            sca_pix = acos(term1 + term2 * term3) * OEL_RADEG;
            binl1c->sca_3D[gd_row][gd_col][view] = binl1c->sca_3D[gd_row][gd_col][view] + sca_pix;
            rotangle = rot_angle(l1rec->senz[pix], l1rec->solz[pix], l1rec->sena[pix],
                                 l1rec->sola[pix]);  // in degrees
            binl1c->rot_angle[gd_row][gd_col][view] = binl1c->rot_angle[gd_row][gd_col][view] + rotangle;
        }

        ibp = pix * nbands;

        for (sb = 0; sb < nbands; sb++) {
           if (view != BAD_FLT && tilt != fill_tilt && gd_row >= 0 && gd_col >= 0 &&
                gd_row <= binl1c->num_gridlines - 1 && gd_col <= binl1c->nbinx - 1 && l1rec->Lt[ibp] != binl1c->fillval2 && l1rec->Lt[ibp] >= minval2_all[0] && 10 * l1rec->Lt[ibp] <= maxval2_all[0]){
                binl1c->I_4D[gd_row][gd_col][view][sb] =
                    binl1c->I_4D[gd_row][gd_col][view][sb] + 10 * l1rec->Lt[ibp];    
                binl1c->nrec_4D_band[gd_row][gd_col][view][sb] += 1;
                binl1c->nrec_3D_view[gd_row][gd_col][view] += 1;

                if(binl1c->i_diff2[gd_row][gd_col][view][sb]==0)
                {  
                     seed_mean=10 * l1rec->Lt[ibp];
                     binl1c->i_diff2[gd_row][gd_col][view][sb]=0.01;//offset substracted before stdev
                }
                else{
                 binl1c->i_diff2[gd_row][gd_col][view][sb]+=(10 * l1rec->Lt[ibp]-seed_mean)*(10 * l1rec->Lt[ibp]-seed_mean);
                }

                if (gd_row > binl1c->num_gridlines - 1 || gd_col > binl1c->nbinx - 1) {
                    if (binl1c->verbose) {
                        cout << "ERROR IN BINNING Lt/rhot --gd_row/gd-col outside the grid limits" << endl;
                        cout << "row.." << gd_row << "gd_col.." << gd_col << endl;
                    }
                }
                // missing QC an QC_bitwise variables
            }
            ibp++;
        }  // bands

    }  // pixels

    return 0;
}

int check_l1c_time(const char *l1b_file, const char *l1c_grid, bin_str *binl1c) {
    string l1c_str = l1c_grid;
    string l1b_str = l1b_file;

    NcFile *nc_l1cgrid;
    try {
        nc_l1cgrid = new NcFile(l1c_grid, NcFile::read);
    } catch (NcException &e) {
        cerr << e.what() << "\n -E- l1cgen l1c_pflag= 8:: Failure reading L1C grid: " + l1c_str << endl;
        exit(EXIT_FAILURE);
    }

    NcDim yd = nc_l1cgrid->getDim("bins_along_track");
    int32_t num_gridlines = yd.getSize();
    double *time_nad_l1c = (double *)calloc(num_gridlines, sizeof(double));
    NcGroup ba_grp = nc_l1cgrid->getGroup("bin_attributes");
    NcVar v1 = ba_grp.getVar("nadir_view_time");
    v1.getVar(time_nad_l1c);

    NcFile *nc_l1b;

    try {
        nc_l1b = new NcFile(l1b_file, NcFile::read);
    } catch (NcException &e) {
        cerr << e.what() << "\n -E- l1cgen l1c_pflag= 8:: Failure reading L1B or L1C-ANCfile: " + l1b_str
             << endl;
        exit(EXIT_FAILURE);
    }

    // only if L1C ANC file
    int32_t nscans;

    if (binl1c->l1c_anc[0] != '\0' && strcmp(binl1c->l1c_anc, l1b_file) == 0) {
        if (binl1c->verbose)
            cout << "L1C ANC file provided" << endl;
        yd = nc_l1b->getDim("bins_along_track");
        int32_t num_gridlines2 = yd.getSize();
        nscans = num_gridlines2;

        if (num_gridlines != num_gridlines2) {
            cerr << "WARNING number of gridlines ARE NOT THE SAME IN L1C grid and ANC files"
                 << "L1C grid # " << num_gridlines << "ANC file # " << num_gridlines2 << endl;
        }
    } else {
        yd = nc_l1b->getDim("scans");
        if(yd.isNull()) {
            yd = nc_l1b->getDim("number_of_scans");
            if(yd.isNull()) {
                cerr << "ERROR - could not read number of scans dimension" << endl;
                exit(EXIT_FAILURE);
            }
        }
        nscans = yd.getSize();
    }

    double tend_l1b, *time_l1b = nullptr;

    if (binl1c->l1c_anc[0] == '\0' || strcmp(binl1c->l1c_anc, l1b_file) != 0) {
        time_l1b = (double *)calloc(nscans, sizeof(double));
        NcGroup sla_grp = nc_l1b->getGroup("scan_line_attributes");
        v1 = sla_grp.getVar("time");
        v1.getVar(time_l1b);
        tend_l1b = time_l1b[nscans - 1];
    } else {
        time_l1b = time_nad_l1c;
        tend_l1b = time_l1b[nscans - 2];
    }

    // check time limits
    double tini_l1c = time_nad_l1c[0];
    double tend_l1c = time_nad_l1c[num_gridlines - 1];
    double tini_l1b = time_l1b[0];

    binl1c->tini_l1c = tini_l1c;
    binl1c->tend_l1c = tend_l1c;
    binl1c->tini_l1b = tini_l1b;
    binl1c->tend_l1b = tend_l1b;

    int time_flag = -1;
    if (tini_l1b > (tini_l1c - 5 * 60 - 1) && tini_l1b < (tend_l1c + 5 * 60 + 1)) {
        time_flag = 0;
        if (binl1c->verbose) {
            cout << "OK--L1B granule is INSIDE of the L1C grid--" << endl;
            cout << "tini_l1c.." << tini_l1c << "tend_l1c.." << tend_l1c << "tini_l1b.." << tini_l1b
                 << "tend_l1b.." << tend_l1b << endl;
            cout << "tini_l1c-5*60" << tini_l1c - 5 * 60 << "tend_l1c+5*60" << tend_l1c + 5 * 60 << endl;
        }
    } else {
        time_flag = 1;
        if (binl1c->verbose) {
            cout << "ERROR--L1B granule is outside of the L1C grid--" << endl;
            cout << "tini_l1c.." << tini_l1c << "tend_l1c.." << tend_l1c << "tini_l1b.." << tini_l1b
                 << "tend_l1b.." << tend_l1b << endl;
            cout << "tini_l1c-5*60" << tini_l1c - 5 * 60 << "tend_l1c+5*60" << tend_l1c + 5 * 60 << endl;
        }
    }

    if (time_nad_l1c != nullptr)
        free(time_nad_l1c);
    time_nad_l1c = nullptr;

    if (binl1c->l1c_anc[0] == '\0' || strcmp(binl1c->l1c_anc, l1b_file) != 0) {
        if (time_l1b != nullptr)
            free(time_l1b);
        time_l1b = nullptr;
    }

    nc_l1cgrid->close();
    nc_l1b->close();

    return time_flag;
}

int open_l1c(const char *ifile_l1c, size_t *ybins, size_t *xbins, float **lat_gd, float **lon_gd) {
    std::string str;
    std::string ifile_str;
    string gridname, azeast_name;
    std::string fname_out, pathstr, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, swtnum, extstr,
        granstr, timestr, azstr, missionstr, ofilestr;

    //  if(binl1c->verbose)cout<<"Opening L1C
    //  grid........................................................................."<<endl;
    string l1c_str(ifile_l1c);

    NcFile *nc_l1cgrid;

    try {
        nc_l1cgrid = new NcFile(ifile_l1c, NcFile::read);
    } catch (NcException &e) {
        e.what();
        cerr << "l1cgen :: Failure reading L1C grid: " + l1c_str << endl;
        exit(1);
    }
    NcGroup geo_grp = nc_l1cgrid->getGroup("geolocation_data");
    NcVar v1 = geo_grp.getVar("latitude");
    v1.getVar(&lat_gd[0][0]);
    v1 = geo_grp.getVar("longitude");
    v1.getVar(&lon_gd[0][0]);

    nc_l1cgrid->close();

    return 0;
}

int open_l1c_grid(const char *ifile_l1c, bin_str *binl1c, float **lat_gd, float **lon_gd, float **alt_gd) {
    std::string str;
    std::string ifile_str;
    string gridname, azeast_name;
    std::string fname_out, pathstr, senstr, monstr, daystr, yearstr, prodstr, gdstr, swtstr, swtnum, extstr,
        granstr, timestr, azstr, missionstr, ofilestr;

    string l1c_str(ifile_l1c);

    NcFile *nc_l1cgrid;

    try {
        nc_l1cgrid = new NcFile(ifile_l1c, NcFile::read);
    } catch (NcException &e) {
        cerr << e.what() << "\n -E- l1cgen :: Failure reading L1C grid: " + l1c_str << endl;
        exit(1);
    }
    NcGroup geo_grp = nc_l1cgrid->getGroup("geolocation_data");
    NcVar v1 = geo_grp.getVar("latitude");
    v1.getVar(&lat_gd[0][0]);
    v1 = geo_grp.getVar("longitude");
    v1.getVar(&lon_gd[0][0]);
    v1 = geo_grp.getVar("height");
    v1.getVar(&alt_gd[0][0]);

    NcDim ydim = nc_l1cgrid->getDim("bins_along_track");
    NcDim xdim = nc_l1cgrid->getDim("bins_across_track");
    binl1c->num_gridlines = ydim.getSize();
    binl1c->nbinx = xdim.getSize();

    nc_l1cgrid->close();

    return 0;
}


double search_calc_dotprod(bin_str *binl1c, float bvec[3], int32_t gridline) {
    int32_t nbinx = binl1c->nbinx;
    float gnvec[3];
    float gnvm;

    // normal vectors for L1C rows
    if (binl1c->lat_gd[gridline][nbinx - 1] > 90 || binl1c->lat_gd[gridline][nbinx - 1] < -90 ||
        binl1c->lon_gd[gridline][nbinx - 1] < -180 || binl1c->lon_gd[gridline][nbinx - 1] > 180) {
        if (binl1c->verbose)
            cout << "lat lon for L1C GRID out of the boundaries.." << endl;
        exit(1);
    }
    if (binl1c->lat_gd[gridline][0] > 90 || binl1c->lat_gd[gridline][0] < -90 || binl1c->lon_gd[gridline][0] < -180 ||
        binl1c->lon_gd[gridline][0] > 180) {
        if (binl1c->verbose)
            cout << "lat lon for L1C GRID out of the boundaries.." << endl;
        exit(1);
    }

    gnvec[0] = sin(binl1c->lon_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    cos(binl1c->lat_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    sin(binl1c->lat_gd[gridline][0] * OEL_DEGRAD) -
                sin(binl1c->lat_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    sin(binl1c->lon_gd[gridline][0] * OEL_DEGRAD) * cos(binl1c->lat_gd[gridline][0] * OEL_DEGRAD);
    gnvec[1] = sin(binl1c->lat_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    cos(binl1c->lon_gd[gridline][0] * OEL_DEGRAD) * cos(binl1c->lat_gd[gridline][0] * OEL_DEGRAD) -
                cos(binl1c->lon_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    cos(binl1c->lat_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    sin(binl1c->lat_gd[gridline][0] * OEL_DEGRAD);
    gnvec[2] = cos(binl1c->lon_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    cos(binl1c->lat_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    sin(binl1c->lon_gd[gridline][0] * OEL_DEGRAD) * cos(binl1c->lat_gd[gridline][0] * OEL_DEGRAD) -
                sin(binl1c->lon_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    cos(binl1c->lat_gd[gridline][nbinx - 1] * OEL_DEGRAD) *
                    cos(binl1c->lon_gd[gridline][0] * OEL_DEGRAD) * cos(binl1c->lat_gd[gridline][0] * OEL_DEGRAD);

    // vector norm
    gnvm = sqrt(gnvec[0] * gnvec[0] + gnvec[1] * gnvec[1] + gnvec[2] * gnvec[2]);
    if (isnan(gnvm) == 1) {
        if (binl1c->verbose)
            cout << "NAN value for gnvm.." << endl;
        exit(1);
    }
    if (gnvm == 0) {
        if (binl1c->verbose)
            cout << "ERROR gnvm == 0--- WE CANT NORMALIZE..." << endl;
        exit(1);
    }
    // normalization
    gnvec[0] = gnvec[0] / gnvm;
    gnvec[1] = gnvec[1] / gnvm;
    gnvec[2] = gnvec[2] / gnvm;

    // dot prod, orbital normaliz and transposed by pixel vector
    return gnvec[0] * bvec[0] + gnvec[1] * bvec[1] + gnvec[2] * bvec[2];
}




int search_l1c_parallax(filehandle *l1file, float *lat_new, float *lon_new, bin_str *binl1c,
                        short **gdindex) {
    int32_t num_gridlines, nbinx;   
    int flag_out = -1;
    size_t pix,badpix;
    int32_t i;
    size_t num_pixels;
    int irow = -1, col = -1;
    float bmcm, bm = 100;
    double db;
    float c1, c2, c3;
    double dotprod, dot_firstline, dot_lastline;
    int32_t last_dotprod_index = 0;
    double last_dotprod;
    float gvec[3], bvec[3];

    num_gridlines = binl1c->num_gridlines;
    nbinx = binl1c->nbinx;
    num_pixels = l1file->npix; 
    db = (5.2) / 6371 / 2;  // Half of bin size in radians, resolution in km
    double db_fudge = db * 1.5;

    // big loop
    // dot product
    for (pix = 0; pix < num_pixels; pix++) {
        gdindex[pix][0] = -1;  // actually these are row/col not irow/icol like in my searching algos!!!
        gdindex[pix][1] = -1;

        if (lat_new[pix] > 90 || lat_new[pix] < -90 || lon_new[pix] < -180 || lon_new[pix] > 180 ||
            lat_new[pix] == BAD_FLT || lon_new[pix] == BAD_FLT) {
            if (binl1c->verbose)
                cout
                    << "ERROR in search_l1cgen --latitude longitude pixel out of the boundaries.OR FILLVALUE."
                    << "latpix>90 or <-90.." << lat_new[pix] << "lonpix>180 or <-180.." << lon_new[pix]
                    << endl;
            flag_out = 110;
            return flag_out;
        }
        bvec[0] = cos(lon_new[pix] * OEL_DEGRAD) * cos(lat_new[pix] * OEL_DEGRAD);
        bvec[1] = sin(lon_new[pix] * OEL_DEGRAD) * cos(lat_new[pix] * OEL_DEGRAD);
        bvec[2] = sin(lat_new[pix] * OEL_DEGRAD);

        last_dotprod = search_calc_dotprod(binl1c, bvec, last_dotprod_index);
        if(last_dotprod < 0)
            last_dotprod = 0 - last_dotprod;
        double min_dotprod = last_dotprod;
        int min_dotprod_index = last_dotprod_index;
	    bool found_new = false;

        // look right
        for (i = last_dotprod_index+1; i < num_gridlines; i++) {
            dotprod = search_calc_dotprod(binl1c, bvec, i);
            if(dotprod < 0)
                dotprod = 0 - dotprod;
            if(dotprod > last_dotprod)
                break;
            if (dotprod < min_dotprod) {
                found_new = true;
                min_dotprod = dotprod;
                min_dotprod_index = i;
                last_dotprod = dotprod;
                last_dotprod_index = i;
            }
        }  // end lines

        // look left
        if(!found_new && last_dotprod_index > 0) {
            for (i = last_dotprod_index-1; i >= 0; i--) {
                dotprod = search_calc_dotprod(binl1c, bvec, i);
                if(dotprod < 0)
                    dotprod = 0 - dotprod;
                if(dotprod > last_dotprod)
                    break;
                if (dotprod < min_dotprod) {
                    min_dotprod = dotprod;
                    min_dotprod_index = i;
                    last_dotprod = dotprod;
                    last_dotprod_index = i;
                }
            }  // end lines
        } // not found_new

        if (min_dotprod <= db_fudge) {
            gdindex[pix][0] = min_dotprod_index + 1;
        }

        // calc first gridline
        dot_firstline = search_calc_dotprod(binl1c, bvec, 0);

        // calc last gridline
        dot_lastline = search_calc_dotprod(binl1c, bvec, num_gridlines - 1);

        // for each pixels
        if (dot_firstline <= db && dot_lastline > -db) {
            // find column
            irow = gdindex[pix][0] - 1;

            if (irow < 0) {
                if (binl1c->verbose)cout << "ERROR irow<0 in search_l1c..." << irow << "at pix#.." << pix + 1 << endl;
                flag_out = 110;
                return flag_out;
            }

            for (int j = 0; j < nbinx; j++) {
                gvec[0] =
                    cos(binl1c->lon_gd[irow][j] * OEL_DEGRAD) * cos(binl1c->lat_gd[irow][j] * OEL_DEGRAD);
                gvec[1] =
                    sin(binl1c->lon_gd[irow][j] * OEL_DEGRAD) * cos(binl1c->lat_gd[irow][j] * OEL_DEGRAD);
                gvec[2] = sin(binl1c->lat_gd[irow][j] * OEL_DEGRAD);

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
                if (binl1c->verbose)
                    cout << "ERROR col<1 in search_l1c..." << col << "at pix#.." << pix + 1 << endl;
                flag_out = 110;
                return flag_out;
            }
            gdindex[pix][1] = col;

            if (irow == num_gridlines - 1 && lat_new[pix] < 85.0) {
                gdindex[pix][0] = -1;
                gdindex[pix][1] = -1;
            }

            bm = 100;
            col = -1;
        } else {
            gdindex[pix][0] = -1;  // actually these are row/col not irow/icol like in my searching algos!!!
            gdindex[pix][1] = -1;
        }
    }  // end pixels

    for (pix = 0; pix < num_pixels; pix++) {
        if (gdindex[pix][0] < 1 && gdindex[pix][1] < 1) {
            badpix++;
        }
    }

    if (badpix == num_pixels) {
        flag_out = 110;
		 } else
        flag_out = 0;

    if (binl1c->verbose && flag_out == 110)
        cout << "THIS LINE WILL BE SKIPPED -- NO PIXELS BINNED.............." << endl;

   return flag_out;

}

int search_l1c(filehandle *l1file, l1str *l1rec, bin_str *binl1c, short **gdindex) {
    int32_t num_gridlines, nbinx;
    float gvec[3], bvec[3];
    int flag_out = 0;
    size_t pix;
    size_t badpix = 0;
    int32_t i;
    size_t num_pixels;
    int irow = -1, col = -1;
    float bmcm, bm = 100;
    double db;
    float c1, c2, c3;
    double dotprod, dot_firstline, dot_lastline;

    int32_t last_dotprod_index = 0;
    double last_dotprod;

    num_gridlines = binl1c->num_gridlines;
    nbinx = binl1c->nbinx;
    num_pixels = l1file->npix;

    db = (5.2) / 6371 / 2;  // Half of bin size in radians, resolution in km
    double db_fudge = db * 1.5;

    // big loop
    // dot product
    for (pix = 0; pix < num_pixels; pix++) {
        gdindex[pix][0] = -1;  // actually these are row/col not irow/icol like in my searching algos!!!
        gdindex[pix][1] = -1;
        flag_out = 0;

        if (l1rec->lat[pix] > 90 || l1rec->lat[pix] < -90 || l1rec->lon[pix] < -180 ||
            l1rec->lon[pix] > 180 || l1rec->lat[pix] == BAD_FLT || l1rec->lon[pix] == BAD_FLT) {
            if (binl1c->verbose)
                cout
                    << "ERROR in search_l1cgen --latitude longitude pixel out of the boundaries.OR FILLVALUE."
                    << "latpix>90 or <-90.." << l1rec->lat[pix] << "lonpix>180 or <-180.." << l1rec->lon[pix]
                    << endl;
            binl1c->badgeo++;
            flag_out = 110;
            return flag_out;
        }

        bvec[0] = cos(l1rec->lon[pix] * OEL_DEGRAD) * cos(l1rec->lat[pix] * OEL_DEGRAD);
        bvec[1] = sin(l1rec->lon[pix] * OEL_DEGRAD) * cos(l1rec->lat[pix] * OEL_DEGRAD);
        bvec[2] = sin(l1rec->lat[pix] * OEL_DEGRAD);

        last_dotprod = search_calc_dotprod(binl1c, bvec, last_dotprod_index);
        if(last_dotprod < 0)
            last_dotprod = 0 - last_dotprod;
        double min_dotprod = last_dotprod;
        int min_dotprod_index = last_dotprod_index;
        bool found_new = false;

        // look right
        for (i = last_dotprod_index+1; i < num_gridlines; i++) {
            dotprod = search_calc_dotprod(binl1c, bvec, i);
            if(dotprod < 0)
                dotprod = 0 - dotprod;
            if(dotprod > last_dotprod)
                break;
            if (dotprod < min_dotprod) {
                found_new = true;
                min_dotprod = dotprod;
                min_dotprod_index = i;
                last_dotprod = dotprod;
                last_dotprod_index = i;
            }
        }  // end lines

        // look left
        if(!found_new && last_dotprod_index > 0) {
            for (i = last_dotprod_index-1; i >= 0; i--) {
                dotprod = search_calc_dotprod(binl1c, bvec, i);
                if(dotprod < 0)
                    dotprod = 0 - dotprod;
                if(dotprod > last_dotprod)
                    break;
                if (dotprod < min_dotprod) {
                    min_dotprod = dotprod;
                    min_dotprod_index = i;
                    last_dotprod = dotprod;
                    last_dotprod_index = i;
                }
            }  // end lines
        } // not found_new

        if (min_dotprod <= db_fudge) {
            gdindex[pix][0] = min_dotprod_index + 1;
        }

        // calc first gridline
        dot_firstline = search_calc_dotprod(binl1c, bvec, 0);

        // calc last gridline
        dot_lastline = search_calc_dotprod(binl1c, bvec, num_gridlines - 1);

        // for each pixels
        if (dot_firstline <= db && dot_lastline > -db) {
            // find column
            irow = gdindex[pix][0] - 1;

            if (irow < 0) {
                if (binl1c->verbose)
                    cout << "ERROR irow<0 in search_l1c..." << irow << "at pix#.." << pix + 1 << endl;
                flag_out = 110;
                return flag_out;
            }

            for (int j = 0; j < nbinx; j++) {
                gvec[0] =
                    cos(binl1c->lon_gd[irow][j] * OEL_DEGRAD) * cos(binl1c->lat_gd[irow][j] * OEL_DEGRAD);
                gvec[1] =
                    sin(binl1c->lon_gd[irow][j] * OEL_DEGRAD) * cos(binl1c->lat_gd[irow][j] * OEL_DEGRAD);
                gvec[2] = sin(binl1c->lat_gd[irow][j] * OEL_DEGRAD);

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
                if (binl1c->verbose)
                    cout << "ERROR col<1 in search_l1c..." << col << "at pix#.." << pix + 1 << endl;
                flag_out = 110;
                return flag_out;
            }
            gdindex[pix][1] = col;

            if (irow == num_gridlines - 1 && l1rec->lat[pix] < 85.0) {
                gdindex[pix][0] = -1;
                gdindex[pix][1] = -1;
            }

            bm = 100;
            col = -1;
        } else {
            gdindex[pix][0] = -1;  // actually these are row/col not irow/icol like in my searching algos!!!
            gdindex[pix][1] = -1;
        }
    }  // end pixels

    for (pix = 0; pix < num_pixels; pix++) {
        if (gdindex[pix][0] < 1 && gdindex[pix][1] < 1) {
            badpix++;
        }
    }

    if (badpix == num_pixels) {
        flag_out = 110;
    } else
        flag_out = 0;

    if (binl1c->verbose && flag_out == 110)
        cout << "THIS LINE WILL BE SKIPPED -- NO PIXELS BINNED.............." << endl;

    return flag_out;
}

// fred/richard implemtation
// one pixel at the time----
int search2_l1c(size_t ybins, size_t xbins, float lat, float lon, float **lat_gd, float **lon_gd, short *row,
                short *col) {
    int32_t num_gridlines, nbinx, ic;
    short gdrow = -1, gdcol = -1;
    float pvec[3], gnvm, db;
    float **gnvec = nullptr, **cnvec = nullptr, **dc = nullptr;
    int flag_out = -1;
    int32_t i;
    float dotprod, dotprod2, dcm;
    flag_out = -1;

    num_gridlines = ybins;
    nbinx = xbins;

    db = (5.2) / 6371 / 2;  // Half of bin size in radians, resolution in km

    // fred/richard implemtation
    gnvec = allocate2d_float(3, num_gridlines);
    cnvec = allocate2d_float(3, num_gridlines);
    dc = allocate2d_float(3, num_gridlines);

    ic = nbinx / 2;

    for (i = 0; i < num_gridlines; i++) {
        // normal vectors for L1C rows
        if (lat_gd[i][nbinx - 1] > 90 || lat_gd[i][nbinx - 1] < -90 || lon_gd[i][nbinx - 1] < -180 ||
            lon_gd[i][nbinx - 1] > 180) {
            exit(1);
        }
        if (lat_gd[i][0] > 90 || lat_gd[i][0] < -90 || lon_gd[i][0] < -180 || lon_gd[i][0] > 180) {
            exit(1);
        }
        // compute normal vector for each rowgrid----
        gnvec[0][i] = sin(lon_gd[i][nbinx - 1] * OEL_DEGRAD) * cos(lat_gd[i][nbinx - 1] * OEL_DEGRAD) *
                          sin(lat_gd[i][0] * OEL_DEGRAD) -
                      sin(lat_gd[i][nbinx - 1] * OEL_DEGRAD) * sin(lon_gd[i][0] * OEL_DEGRAD) *
                          cos(lat_gd[i][0] * OEL_DEGRAD);
        gnvec[1][i] = sin(lat_gd[i][nbinx - 1] * OEL_DEGRAD) * cos(lon_gd[i][0] * OEL_DEGRAD) *
                          cos(lat_gd[i][0] * OEL_DEGRAD) -
                      cos(lon_gd[i][nbinx - 1] * OEL_DEGRAD) * cos(lat_gd[i][nbinx - 1] * OEL_DEGRAD) *
                          sin(lat_gd[i][0] * OEL_DEGRAD);
        gnvec[2][i] = cos(lon_gd[i][nbinx - 1] * OEL_DEGRAD) * cos(lat_gd[i][nbinx - 1] * OEL_DEGRAD) *
                          sin(lon_gd[i][0] * OEL_DEGRAD) * cos(lat_gd[i][0] * OEL_DEGRAD) -
                      sin(lon_gd[i][nbinx - 1] * OEL_DEGRAD) * cos(lat_gd[i][nbinx - 1] * OEL_DEGRAD) *
                          cos(lon_gd[i][0] * OEL_DEGRAD) * cos(lat_gd[i][0] * OEL_DEGRAD);

        gnvm = sqrt(gnvec[0][i] * gnvec[0][i] + gnvec[1][i] * gnvec[1][i] + gnvec[2][i] * gnvec[2][i]);
        if (isnan(gnvm) == 1 || gnvm == 0) {
            exit(1);
        }
        // normalization
        gnvec[0][i] = gnvec[0][i] / gnvm;
        gnvec[1][i] = gnvec[1][i] / gnvm;
        gnvec[2][i] = gnvec[2][i] / gnvm;

        // Compute normals to center columns
        cnvec[0][i] = sin(lon_gd[i][ic] * OEL_DEGRAD) * cos(lat_gd[i][ic] * OEL_DEGRAD) * gnvec[2][i] -
                      sin(lat_gd[i][ic] * OEL_DEGRAD) * gnvec[1][i];
        cnvec[1][i] = sin(lat_gd[i][ic] * OEL_DEGRAD) * gnvec[0][i] -
                      cos(lon_gd[i][ic] * OEL_DEGRAD) * cos(lat_gd[i][ic] * OEL_DEGRAD) * gnvec[2][i];
        cnvec[2][i] = cos(lon_gd[i][ic] * OEL_DEGRAD) * cos(lat_gd[i][ic] * OEL_DEGRAD) * gnvec[1][i] -
                      sin(lon_gd[i][ic] * OEL_DEGRAD) * cos(lat_gd[i][ic] * OEL_DEGRAD) * gnvec[0][i];

        //; Compute grid row nadir resolution
        dc[0][i] = cos(lon_gd[i][ic + 1] * OEL_DEGRAD) * cos(lat_gd[i][ic + 1] * OEL_DEGRAD) -
                   cos(lon_gd[i][ic] * OEL_DEGRAD) * cos(lat_gd[i][ic] * OEL_DEGRAD);
        dc[1][i] = sin(lon_gd[i][ic + 1] * OEL_DEGRAD) * cos(lat_gd[i][ic + 1] * OEL_DEGRAD) -
                   sin(lon_gd[i][ic] * OEL_DEGRAD) * cos(lat_gd[i][ic] * OEL_DEGRAD);
        dc[2][i] = sin(lat_gd[i][ic + 1] * OEL_DEGRAD) - sin(lat_gd[i][ic] * OEL_DEGRAD);
    }

    // dot product
    for (i = 0; i < num_gridlines; i++) {
        if (lat > 90 || lat < -90 || lon < -180 || lon > 180) {
            exit(1);
        }

        pvec[0] = cos(lon * OEL_DEGRAD) * cos(lat * OEL_DEGRAD);
        pvec[1] = sin(lon * OEL_DEGRAD) * cos(lat * OEL_DEGRAD);
        pvec[2] = sin(lat * OEL_DEGRAD);

        dotprod = pvec[0] * gnvec[0][i] + pvec[1] * gnvec[1][i] + pvec[2] * gnvec[2][i];

        // check row for pixel j
        if (dotprod <= db && dotprod > -db && gdrow < 0) {
            gdrow = i;
        }

        if (gdrow >= 0) {
            dotprod2 = pvec[0] * cnvec[0][gdrow] + pvec[1] * cnvec[1][gdrow] + pvec[2] * cnvec[2][gdrow];
            dcm =
                sqrt(dc[0][gdrow] * dc[0][gdrow] + dc[1][gdrow] * dc[1][gdrow] + dc[2][gdrow] * dc[2][gdrow]);
            gdcol = ic + dotprod2 / dcm + .5;
            if (gdrow >= 0 && gdcol >= 0) {
                flag_out = 0;
            } else {
                flag_out = 1;
            }
        }
    }
    *row = gdrow;
    *col = gdcol;

    gdrow = -1;
    gdcol = -1;

    // bmcm only checked one pixel!!I
    if (dc != nullptr)
        free2d_float(dc);
    if (gnvec != nullptr)
        free2d_float(gnvec);
    if (cnvec != nullptr)
        free2d_float(cnvec);

    if (flag_out == 0)
        return 0;
    else
        return 1;
}

double rot_angle(double senz, double solz, double sena, double suna) {
    double rotangle = double(BAD_FLT), term2;
    if (senz != BAD_FLT && solz != BAD_FLT && sena != BAD_FLT && suna != BAD_FLT) {
        double cosunz = cos(solz * OEL_DEGRAD);
        double cosenz = cos(senz * OEL_DEGRAD);
        double term1 = -1 * cosunz + cosenz * cos(senz * OEL_DEGRAD + OEL_DEGRAD);
        double cos2 = 0.5 * (cos(2 * senz * OEL_DEGRAD) + 1);
        double term3 =
            sqrt(1 - cos(senz * OEL_DEGRAD + OEL_DEGRAD) * cos(senz * OEL_DEGRAD + OEL_DEGRAD));

        if ((sena - suna) < 2 * 180. && (sena - suna) > 180.)
            term2 = sqrt(1 - cos2) * term3;
        if ((sena - suna) < 180. && (sena - suna) > 0.)
            term2 = -1 * sqrt(1 - cos2) * term3;

        double cos_alpha = term1 / term2;

        rotangle = acos(cos_alpha) * OEL_RADEG;  // in degrees
    }

    return rotangle;
}
