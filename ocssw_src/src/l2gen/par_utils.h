#ifndef _PARUTILS_H
#define _PARUTILS_H

#include <math.h>
#include <stdio.h>

#include "filehandle.h"
#include "l12_proto.h"
#include "l2_struc.h"

// typedef int logical;
//  #define FALSE 0
//  #define TRUE 1
#define NMODEL 12
#define NPHASE 75

#define EARTH_SURF_PRESSURE 1013.15f
#define TAU_550_LOW_BOUND 0.1f
#define ANGSTROM_DEFAULT_VALUE 0.2537f
#define MAXGLINT 0.05
#define MAXHOURS 48
#define EINSTEIN 1.193
#define MAX_SOLZEN 89.5f
#define PAR_0 176.323f
#define PAR_0_2023 176.41f
#define TAUCONS_LOW 15.0f
#define TAUCONS_HIGH 300.0f

/**
 * @brief Binary search algorithm
 *
 * @param arr input array, sorted (either desending or ascending order. It is
 * user responsibility to sort the array)
 * @param s start
 * @param e end
 * @param val value
 * @return int index to find
 */
size_t search(const float *arr, size_t s, size_t e, float val, size_t *i_s,
              size_t *i_e);

/**
 * @brief Kasten eq to compute airmass
 *
 * @param solz - angle
 * @return airmass
 */
float kasten_equation(float solz);

typedef struct Temp_merra_data {
    float *wind, *tau550, *angstrom, *surfPress, *dobson, *watVap, *tauCld, *cf,
            *lat, *lon;
    float start_time, end_time, obs_time, time_rise_min, time_set_max;
    size_t tstep, dim_x, dim_y;
} merra2_temp;

/**
 * @brief N-dimensional interpolation
 *
 * @param n_dims - number of dimensions
 * @param dims - arrays which contains dimension sizes
 * @param point - point of interest ( array of size n_dims )
 * @param grid - two dimensional array, where grid[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @param lut - look up table, where  lut[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @return float
 */
float interpnd(size_t n_dims, const size_t *dims, const float *point,
               float **grid, const float *lut);

/**
 * @brief 4-dimensional interpolation
 *
 *  @param dims - arrays which contains dimension sizes
 * @param point - point of interest ( array of size n_dims )
 * @param grid - two dimensional array, where grid[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @param lut - look up table, where  lut[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @return float
 */
float interp4d(const size_t *dims, const float *point, float **grid,
               const float *lut);

/**
 * @brief 3-dimensional interpolation
 *
 *  @param dims - arrays which contains dimension sizes
 * @param point - point of interest ( array of size n_dims )
 * @param grid - two dimensional array, where grid[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @param lut - look up table, where  lut[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @return float
 */
float interp3d(const size_t *dims, const float *point, float **grid,
               const float *lut);

/**
 * @brief 1-dimensional interpolation
 *
 *  @param dims - arrays which contains dimension sizes
 * @param point - point of interest ( array of size n_dims )
 * @param grid - two dimensional array, where grid[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @param lut - look up table, where  lut[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @return float
 */
float interp1d(const size_t *dims, const float *point, float **grid,
               const float *lut);

/**
 * @brief 6-dimensional interpolation
 *
 *  @param dims - arrays which contains dimension sizes
 * @param point - point of interest ( array of size n_dims )
 * @param grid - two dimensional array, where grid[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @param lut - look up table, where  lut[i_dim] is an array of  size
 * dims[i_dim] where 0<= i_dim < n_dims
 * @return float
 */
float interp6d(const size_t *dims, const float *point, float **grid,
               const float *lut);

/**
 * @brief look up tables for atmosphere intrinsic reflectance (rho), total
 * gaseous absorption (tg) and total scattering transmittance (td)
 *
 */
typedef struct Rho_dims {
    float *solar_zenith_angle;
    float *view_zenith_angle;
    float *relative_azimuth_angle;
    float *angstrom_coefficients;
    float *optical_depth_at_550_nm;
    float *wavelength;

    size_t dimsolar_zenith_angle;
    size_t dimview_zenith_angle;
    size_t dimrelative_azimuth_angle;
    size_t dimangstrom_coefficients;
    size_t dimoptical_depth_at_550_nm;
    size_t dimwavelength;
} rho_dims;

typedef struct Td_dims {
    float *wavelength;
    float *air_mass;
    float *angstrom_coefficients;
    float *optical_depth_at_550_nm;

    size_t dimwavelength;
    size_t dimair_mass;
    size_t dimangstrom_coefficients;
    size_t dimoptical_depth_at_550_nm;
} td_dims;

typedef struct Tg_dims {
    float *wavelength;
    float *air_mass;
    float *ozone_concentration;
    float *water_vapor_pressure;

    size_t dimwavelength;
    size_t dimair_mass;
    size_t dimozone_concentration;
    size_t dimwater_vapor_pressure;
} tg_dims;

typedef struct Dobson_dims {
    float *days;
    float *latitude;
    size_t dimdays;
    size_t dimlatitude;
} dobson_dims;

typedef struct WaterVapor_dims {
    float *days;
    float *latitude;
    size_t dimdays;
    size_t dimlatitude;
} watervap_dims;

typedef struct Surface_Ocean_Alebedo {
    float *jwl, *achl, *wv1, *bw1, *wv3, *aw3, *wlt_reft, *reft;
    size_t dim_jwvl, dim_wv1, dim_wv3, dim_wlt_reft;
} soa_luts;

typedef struct Scalar_PAR_MU_Cosine {
    float *wind_speed, *doy, *latitude, *CF, *mu_cosine_c, *mu_cosine, *PARo,
            *PARo_c,*PARc_above;
    size_t dim_wind_speed, dim_doy, dim_latitude;

} scalar_par_luts;

typedef struct Scalar_PAR_Inst {
    float *wind_speed, *cos_solz, *cot, *cf_pard_p, *cf_pard_m, *pard_p_cs, *pard_m_cs, *pard_m_oc;
    size_t dim_wind_speed, dim_solz, dim_cot;
} scalar_inst_par_luts;

typedef struct LUTs_data {
    float *lut_rho;
    float *lut_tg;
    float *lut_td;
    rho_dims rhodims;
    td_dims tddims;
    tg_dims tgdims;
    // small LUTs
    float *lut_dobson;
    dobson_dims ozonedims;
    float *lut_watvap;
    watervap_dims watvapdims;
    soa_luts soa_lut;
    scalar_par_luts scalar_luts;
    scalar_inst_par_luts scalar_inst_luts;
} luts_par;

/**
 * @brief
 * //
##############################################################################
// # from Fitzpatrick et al. (2005), compute cloud albedo
// # Fitzpatrick, M. F., R. E. Brandt, and S. G. Warren, 2004: Transmission of
// Solar # Radiation by Clouds over Snow and Ice Surfaces: A Parameterization in
// Terms of # Optical Depth, Solar Zenith Angle, and Surface Albedo. Journal of
// Climate 17, 2: # 266-275,
// https://doi.org/10.1175/1520-0442(2004)017<0266:TOSRBC>2.0.CO;2
//
##############################################################################
 * @param TauCld - MERRA 2 input,cloud optical thickness for a given pixel, a 1D
array resolved with respect to time, hourly
 * @param CF - MERRA 2 input,cloud fraction for a given pixel, aa 1D resolved
with respect to time, hourly
 * @param cosSZ - cosine of solar zenith angle
 * @param t_obs - observed time
 * @param t_range - a 1D time array, the same length as TauCld and CF
 * @param albe_obs - observed albedo for each wavelentgh
 * @param TauCld_obs - observed optical thinckness
 * @param CF_obs - obesrved cloud fraction
 * @param t_step - size of t_range, TauCld and CF
 * @param wl - input wavelength, 1D arrays
 * @param bands_num lengh of wl
 *
 */
void getcldalbe(float *TauCld, float *CF, float cosSZ, float t_obs,
                float *t_range, float *albe_obs, float *TauCld_obs,
                float *CF_obs, size_t t_step, float *wl, size_t bands_num);

float getosa(float wl, float sza, float wind, float chl, float fr,
             const luts_par *luts_data);

void get_luts_data(l2str *l2rec, luts_par *luts_data);

float calc_par(l2str *l2rec, int ip, int nbands, float *Lt, float taua,
               float angstrom, float *wl, float *fo, float *ko3,
               float *taumolbar);

void calc_scalar_inst_par(l2str *l2rec, int ip, float par_above_ins,float * par_scalar_ins);

void calc_scalar_par_mean_cosine(l2str *l2rec, int ip, float par_above,
                                 float par_c, float *scalar_par,
                                 float *mean_cosine);

float calc_par_impl_of_2023(l2str *l2rec, int ip, int nbands, float *Lt,
                            float taua, float angstrom, float *wl, float *fo,
                            float *ko3, float *taumolbar, float *parb,
                            float *parc);

void GetAerPhase(l2str *l2rec, int ip, int32_t nbands, float angstrom,
                 float *phasea, float *omegaa, float *modelAngstrom);

void read_aerosol_par(l2str *l2rec, int32_t nbands, float *tablewavelengths,
                      float *tablephaseangles, float *tablealphas,
                      float *tableomegas, float *tableaerphasefunc);

void *allocateMemoryPar(size_t numBytes, const char *name);

float EstimateDobson(int32_t year, int32_t month, int32_t day, float lat);

float EstimateWatVap(int32_t year, int32_t month, int32_t day, float lat);

float varsol(int32_t jday);

void triseset(int32_t jday, float xlon, float xlat, float *trise, float *tset);

int Greg2Jul(int32_t year, int32_t month, int32_t day);

// return solar zenith - ignore solar azimuth angle
float get_solz(int jday, float time, float lon, float lat);

float interp_as_taulow(float csz, float tau);

float interp_as_tauhigh(float csz);

/**
 * @brief
########################################################################
## Compute sun glint reflectance using wind speed only (no direction)
#######################################################################
 *
 * @param sz -solar zenith angle
 * @param vz - view zenith angle
 * @param ra - radiance
 * @param ws - windspeed
 * @return reflectance glint
 */
float SunGlint(float sz, float vz, float ra, float ws);

#endif
