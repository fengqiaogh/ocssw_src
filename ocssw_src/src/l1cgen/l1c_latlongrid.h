#ifndef L1C_LATLONGRID_H
#define L1C_LATLONGRID_H

#include "filehandle.h"
#include "l1.h"
#include <netcdf.h>
#include <unistd.h>
#include <nc4utils.h>
#include <netcdf>

#ifdef __cplusplus
extern "C" {
#endif

class bin_str {
   protected:
   public:
    bin_str();
    virtual ~bin_str();
    virtual int alloc_bin(bin_str *binstr);
    virtual int close_bin(bin_str *binstr);

    // global attributes
    std::string version;
    std::string history;
    std::string pversion;
    std::string doi;
    int sensor;

    // dims
    int16_t num_gridlines;
    int16_t nbinx;
    int16_t nviews;
    int16_t nbands;
    int32_t num_pixels;
    int32_t nscans;
    bool verbose;

    // geo
    int gdshift;
    float **lat_gd;
    float **lon_gd;
    double ***time_l1b;
    double *time_gd;
    float **alt;

    // binned vars
    short **nrec_2D;   // row/col
    short ***nrec_3D;  // row/col/view
    short ***nrec_3D_view;
    float ****nrec_4D_band;  // row/col/view/bands
    float **alt_mean;
    float **alt_rmse;
    short **alt_2D;
    short **alt_diff2;  // row/col
    float ***suna_3D;
    float ***sunz_3D;
    float ***sena_3D;
    float ***senz_3D;
    float ***sca_3D;  // row/col/view
    float ***rot_angle;
    float ****QC_bitwise_4D;  // row/col/view/bands
    float ***QC_4D;           // row/col/view
    float ****I_4D;           // row/col/view/bands
    float ****i_diff2; 
    float ****I_noise_4D;     // #pixels

    // OCIS line by line
    short *obs_per_view;
    float *QC_bitwise;
    float *QC;
    float *I;
    float *I_noise;

    file_format format;

    filehandle *l1file;
    std::string full_l1cgrid;
    std::string l1cgrid;

    short fillval1 = -32767;
    float fillval2 = -32767.;
    double fillval3 = -32767.;

    int32_t inpix;   // in l1c grid
    int32_t outpix;  // out l1c grid
    int32_t badgeo;  // fillvalue or no geolocation or bad

    double tini_l1c;
    double tend_l1c;
    double tini_l1b;
    double tend_l1b;
    std::string date_mid_grid;

    int bintype;  // 0: discrete, 1: area-weighting

    char outlist[FILENAME_MAX];
    char l1c_anc[FILENAME_MAX];
    int cloudem_flag;
    int cloud_type;  // 0 water and 1 ice
    int dem_flag;//geoid or orthometric =  dem height flag
};

int close_bin(bin_str *binstr);
int alloc_bin(bin_str *binstr);

int meta_l1c_grid(char *gridname, bin_str *binstr, int16_t num_gridlines, netCDF::NcFile *nc_output);
int meta_l1c_global(char *gridname, bin_str *binl1c,int16_t num_gridlines, netCDF::NcFile *nc_output);
int meta_l1c_full(filehandle *l1file, bin_str *binstr, const char *l1c_grid, netCDF::NcFile *nc_output);
int meta_l1c_bin(filehandle *l1file, bin_str *binstr, netCDF::NcFile *nc_output);
int meta_l1c_altvar(bin_str *binstr, netCDF::NcFile *nc_output);  // computes rmse for height

int open_l1c(const char *l1c_grid, size_t *ybins, size_t *xbins, float **lat_gd, float **lon_gd);
int open_l1c_grid(const char *l1c_grid, bin_str *binstr, float **lat_gd, float **lon_gd, float **height_gd);
int search_l1c(filehandle *l1file, l1str *l1rec, bin_str *binstr, short **gdindex);
int search2_l1c(size_t ybins, size_t xbins, float lat, float lon, float **lat_gd, float **lon_gd, short *row,
                short *col);
int search_l1c_parallax(filehandle *l1file, float *lat_new, float *lon_new, bin_str *binl1c, short **gdindex);
int bin_l1c(filehandle *l1file, l1str *l1rec, bin_str *binstr, short **gdindex, netCDF::NcFile *nc_output);
int bintime_l1c(filehandle *l1file, l1str *l1rec, bin_str *binstr, short **gdindex, double scantime,
                netCDF::NcFile *nc_output);
int check_l1c_time(const char *l1b_file, const char *l1c_grid,
                   bin_str *binstr);  // check if chosen l1c grid is correct for the chosen l1b granule used
                                      // for binning and FULL l1c

int rmse_l1c_alt(filehandle *l1file, bin_str *binstr, l1str *l1rec,
                 short **gdindex);  // computes diff squared for rmse for height
double rot_angle(double senz, double solz, double sena, double suna);
int parallax(filehandle *l1file, const char *l1c_anc, const char *l1c_grid, l1str *l1rec, bin_str *binl1c,
             short **gdindex, netCDF::NcFile *nc_output,int32_t sline,int firstcall);
int interp2d_l1c(netCDF::NcFile *nc_output);

#ifdef __cplusplus
}  // extern "C"
#endif

#endif
