/*
  get_ctht.h  will hold the structures for cloud top height storage
*/
#ifndef __CTHT__
#define __CTHT__

#include "l12_proto.h"
#include "l2prod.h"
#include "l1.h"
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>
#include <gsl/gsl_cblas.h>

/*
 *  these keep the lines of cloud heght products, the height, pressure, temp
 *  need to be for a set of lines so that the cloud code can use it
 */
typedef struct ctht_lin_struct
  {
  int32_t npix;
  int32_t iscan;
  //   final product arrays all size npix
  float *cth, *ctp, *ctt;  // the variables needed by cloud microphys props cmp
  unsigned char *ct_phase;     //  phase determined - 0 - water, 1 ice, or 
                               //  BAD_UBYTE for un-retrieved

  //  product for all phases, make band-sequential to simplify output (thats
  //  the way the 3-d products are - band chaanging fastest
  //  only the cost and acost do not have _all to denote all phases
  float *cth_all, *ctp_all, *ctt_all;
  float *cth_cod_all, *cth_lcod_all;
  float *cth_raw, *cth_alb_all, *cost_ss, *acost;
  float *oe_akdiag, *dlcod, *dcod, *dcth, *dctp;
  float *dctt, *dalb;
  int32_t *nitr;
  } ctht_lin_str;

typedef struct ctht_lins_struct
  {
  int32_t nrec;
  int32_t cscan;
  ctht_lin_str **ct;  //  array of structures for each line
  } ctht_lins_str;
/*
 *   uncertainty data, Also known as the Sy
 */
typedef struct ctht_unc_struct
  {
  int32_t nbands, nsurf, nphase; // # bands, land sfc type, cld top phase
  float *icpts;  // Uncertainty fit intercepts
  float *grads;  // Uncertainty fit gradients
  float *fit_rms;  // Uncertainty fit RMS
  float *err_corrs_category;  // Spectral correlation in forward 
                              // model uncertainty
  float noise_err, cal_err;  // noise and calibration error, goes with 
                             // the rest here
  } ctht_unc_str;
/*
 *   LUT data, also called lut, one for each cloud phase
 */
typedef struct ctht_lut_struct
  {
  //  axis sizes
  int32_t lut_nsolz;   // view angle
  int32_t lut_nsenz;
  int32_t lut_nrelaz; 
  int32_t lut_nlev;  // levels
  int32_t lut_npress;  //pressure
  int32_t lut_nalb;  // albedo
  int32_t lut_ncod;  // cloud optical depth
  int32_t lut_naod;  // aerosol optical depth
  int32_t lut_nwave;  // bands
  int32_t lut_ncth;  //  cloud top height
  //  axis values
  float *lut_solz;      // view angle 
  float *lut_senz;
  float *lut_relaz;
                     // levels - none as separate
  float *lut_press;  //pressure
  float *lut_alb;  // albedo
  float *lut_cod;  // cloud optical depth
  float *lut_lcod;  // log cloud optical depth, why not?
  float *lut_aod;  // aerosol optical depth
  float *lut_wave;  // band wavelengths
  float *lut_cth;  //  cloud top height
  //  basic arrays
  float *z_and_pp0;  //  relation between height and relative p - an array
                     //  of size [# lev, 2] with 0 = height, 1 = pressure
  float *lut_pp0;  //  actually derived from z_and_pp0 for each phase table
  float *lut;  //  the spectra f( geom, p, band, cod, cth, sfc albedo )
               //  Lut is sized 
               //  [ solz, senz, relaz, press, band, cod, cth, alb ]
               //  with solz fastest

  // There are many more for the lut, some of questionable need, so later
  } ctht_lut_str;

//  the oe_info_str contains all the inversion results from ams_oe_inversion
//  (however, much more manipulation is required to get ctht)
typedef struct oe_info_struct
  {
  float *x_prd_state;  // the final state vector or the product set found
                       //  this is the 'x'
  gsl_matrix *gm_sx;  //  the sx
  int32_t conv_flag;  // conv
  double cost;  //  j, or j_cost
  double acost;  // a-priori cost or ja
  gsl_matrix *gm_ak;  //  gm_av_kernel or avkernel
  gsl_matrix *gm_gain;  //  gm_gain or G in IDL code
  int32_t nitr;  //  nitr
  
   //  There are more, but we'll make the above be a place holder
  } oe_info_str;
/*
 * define the routines
 */
float *get_ctht_lin( l2str *, int );
int32_t comp_ctht( l2str * );
int32_t comp_ctht_lin( l1str *, ctht_lin_str * );
int32_t init_ctht_parms( ctht_lin_str *, int32_t );
int32_t comp_ctht_tst(l1str *, int32_t, ctht_lin_str * );
int32_t ctht_tbl_init( int32_t, l1str *, ctht_unc_str *, ctht_lut_str *, ctht_lut_str * );
// float ctht_glint( solz, senz, relaz, windsp )
float ctht_glint( float, float, float, float );
int32_t int_4d( float *sxax[4], int32_t snax[4], int64_t n_block, float x[4],
  int64_t[16], float[16], float[4] );
int32_t iint_3d( int32_t[3], double[3], int64_t[9], float[9], float[3] );
gsl_matrix *invert_a_matrix(gsl_matrix *);
// ams_oe_inversion( rhot, gm_sy, xa, gm_sa, tmp_lut, min_cost_loc )
int32_t ams_oe_inversion( double *, gsl_matrix *, double *, gsl_matrix *, 
  double *, int32_t *, int32_t *, oe_info_str * );
int32_t my_lut_funct_3d_oe( gsl_vector *, double *, int32_t *, int32_t,
  double *, double * );
float axis_interp( float *, int32_t, float );
//  these 2 are for some diagnostics
void print_mat_contents(gsl_matrix *matrix, int nrow, int ncol );

#endif
