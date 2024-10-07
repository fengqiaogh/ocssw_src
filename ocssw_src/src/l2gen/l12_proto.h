#ifndef _L12_PROTO_H
#define _L12_PROTO_H

#include <stdint.h>
#include <sys/types.h>
#include <unistd.h>
#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <libgen.h>
#include <string.h>
#include <ctype.h>
#include <netcdf.h>

#include <hdf.h>
#include <mfhdf.h>

#include <l1.h>
#include "input_struc.h"
#include "l12_parms.h"
#include "l1q_struc.h"
#include "l2_struc.h"
#include "target_struc.h"
#include "vcal_struc.h"
#include "filter.h"
#include <timeutils.h>
#include "genutils.h"
#include "aer_struc.h"
#include "l2prod.h"
#include "l2prod_struc.h"
#include "dfutils.h"
#include "allocate2d.h"
#include "epr_api.h"
#include "l2_hdf_generic.h"
#include "flags_iop.h"
#include "table_io_wrapper.h"
#include <clo.h>
#include <uncertainty.h>
#include "navigation.h"
#include "read_l3bin.h"
#include "par_utils.h"
#include "get_nitrate.h"
#include "get_ctht.h"
#include "get_ndvi.h"

#ifdef __cplusplus
extern "C" {
#endif

extern instr *input; /* input parameters structure         */

int getl1rec(int32_t recnum, int32_t dscan, l1str *l1rec);

int loadl1(filehandle *l1file, l1str *l1rec);

int openl2(filehandle *l2file);
int writel2(filehandle *l2file, int32_t recnum, l2str *l2rec, int outfile_number);
int closel2(filehandle *l2file);
//VOIDP scale_sds(float *data, l2prodstr *p, int nbands);
void update_flag_cnts(int32_t *flag_cnt, int32_t* flags, int32_t nflags, int32_t npix, uint32_t init_mask);
void update_flag_cnts16(int32_t *flag_cnt, int16_t* flags, int32_t nflags, int32_t npix, uint32_t init_mask);
void update_qual_cnts(int32_t *flag_cnt, int8_t* flags, int32_t nflags, int32_t npix);
int write_flag_pcnts(idDS ds_id, FILE *fpmeta, int32_t *flag_cnt, int32_t nflags, const char * const flag_lname[], int32_t numScans, int32_t numPixels);
int write_qual_flag_pcnts(idDS ds_id, FILE *fpmeta, int32_t *flag_cnt, int32_t nflags, const char * const flag_lname[]);

int open_target(filehandle *file);
int read_target(filehandle *file, int32_t recnum, tgstr *tgrec);
void close_target(void);

int open_aer(filehandle *file);
int read_aer(filehandle *file, int32_t recnum, aestr *aerec);
void close_aer(filehandle *file);

void init_l2(l2str *l2rec, int32_t nbands);

void free_l1q(void);
int32_t alloc_l2(l1str *l1rec, l2str *l2rec);
void free_l2(l2str *l2rec);
int32_t alloc_target(int32_t npix, int32_t nbands, tgstr *tgrec);
int32_t alloc_aer(int32_t npix, int32_t nbands, aestr *aerec);

void init_l2prod();
l2prodstr *get_l2prod_index(const char *name, int32_t sensorID, int32_t numBands, int32_t numPixels,
        int32_t numScans, int32_t *wave);
void write_product_XML_file(char* filename);

int32_t prodlist(int32_t sensorID, int32_t evalmask, const char *inprod, const char *defprod, char outprod[L1_MAXPROD][32]);

int convl12(l1str *l1rec, l2str *l2rec, int32_t spix, int32_t epix, aestr *aerec);
int convl21(l2str *l2rec, tgstr *tgrec, int32_t spix, int32_t epix, float *Lt, vcstr *vrec);

int32_t l2_seawifs(filehandle *l1file, filehandle *l2file);
int32_t l1b_seawifs(filehandle *l1file, filehandle *ofile,
        int32_t sscan, int32_t escan, int32_t dscan,
        int32_t spixl, int32_t epixl, int32_t dpixl);
int32_t get_modis_calfile(int32_t sd_id, char *file);

void setflagbits_l2(l2str *l2rec, int32_t ipix);  
int setanc(l1str *l1rec);
int acq_sfc_albedo(l1str *l1rec);
int acq_cth_albedo(l1str *l1rec);
int init_cld_dat( l1str *l1rec );
int read_albedo( int32_t d_np, int32_t d_nl, int32_t st_lin, int32_t ix_clim,
  int ncid, int *d_id, unsigned char ***alb_dat );

void cpl1rec(l1str *l1new, l1str *l1old);

int b128_msk_init(char *landfile, int msknum);
int b128_msk_get(float lat, float lon, int msknum);
int dem_init();
int land_mask_init();
int land_mask(float lat, float lon);
int land_bath_mask(l1str *l1rec,int32_t ip);
int bath_mask_init(char *file);
int bath_mask(float lat, float lon);
float get_dem(float lat, float lon);
int get_height(l1str *l1rec, int32_t ip, int terrain_corrected);
void free_deminfo();

void lowercase(char *s);

void atmocor1(l1str *l1rec, int32_t ip);
int atmocor2(l2str *l2rec, aestr *aerec, int32_t ip);
void whitecaps(int32_t sensorID, int32_t evalmask, int32_t nwave, float ws, float ws_max, float rhof[]);
void rayleigh(l1str *l1rec, int32_t ip);
int aerosol(l2str *l2rec, int32_t aer_opt_in, aestr *aerec, int32_t ip,
        float wave[], int32_t nwave, int32_t nir_s_in, int32_t nir_l_in,
        float F0_in[], float La1_in[], float La2_out[],
        float t_sol_out[], float t_sen_out[], float *eps, float taua_out[],
        int32_t *modmin, int32_t *modmax, float *modrat,
        int32_t *modmin2, int32_t *modmax2, float *modrat2, float *mbac_w);
void gaseous_transmittance(l1str *l1rec, int32_t ip);

float ky_airmass(float theta);
float pp_airmass(float theta);

void get_rhown_nir(char *fqfile, float Rrs[], float wave[], int32_t nir_s, int32_t nir_l,
        float aw[], float bbw[], float chl, float rhown[],float dRrs[],float dchl, float drhown[]);
void get_rhown_eval(char *fqfile, float Rrs[], float wave[], int32_t nir_s, int32_t nir_l,
        int32_t nwave, float aw[], float bbw[], float chl,
        float solz, float senz, float phi, float rhown[],l2str *l2rec,int32_t ip);
void get_rho_mumm(l2str *l2rec, int32_t ip, int32_t iw, float *rhom);
void get_rhown_mumm(l2str *l2rec, int32_t ip, int32_t nir_s, int32_t nir_l, float rhown[]);
void glint_rad(int32_t num_iter, int32_t nband, int32_t nir_s, int32_t nir_l,
        float glint_coef, float air_mass,
        float mu0, float F0[], float taur[], float taua[], float La[], float TLg[], uncertainty_t *errstr);

float fresnel_coef(float mu, float n);
float fresnel_sen(float senz, int return_tf);
void fresnel_sol(float wave[], int32_t nwave, float solz, float ws, float brdf[], int return_tf);
void foqint_morel(char *file, float wave[], int32_t nwave, float solz, float senzp,
        float phi, float chl, float brdf[]);
void qint_morel(float wave[], int32_t nwave, float solz, float chl, float Qn[]);
void gothic_R(float wave[], int32_t nwave, float solz, float senz, float ws, float R[]);
int ocbrdf(l2str *l2rec, int32_t ip, int32_t brdf_opt, float wave[], int32_t nwave,
        float solz, float senz, float phi, float ws, float chl, float nLw[], float Fo[], float brdf[]);

void nlw_outband(int32_t evalmask, int32_t sensorID, float wave[], int32_t nwave, float Lw[], float nLw[], float outband_correction[]);

int l2gen_usage(const char *prog);
int msl12_input_defaults(filehandle *l1file);
void msl12_input_init();
int l2gen_init_options(clo_optionList_t* list, const char* prog);
int msl12_option_input(int argc, char *argv[], clo_optionList_t* list,
        char *progName, filehandle *l1file);
int msl12_input(int argc, char *argv[], const char* progName, filehandle *l1file);

  //int windex_get(int32_t wave);

float bin_climatology(char *l3file, int32_t day, float lon, float lat, char *prodname);

float get_default_chl(l2str *l2rec, float Rrs[]);
void get_chl(l2str *l2rec, int prodnum, float prod[]);
void get_las(l2str *l2rec, l2prodstr *p, float prod[]);
//void get_sma(l2str *l2rec, int32_t prodID, float prod[]);
void get_tsm(l2str *l2rec, int prodnum, float prod[]);
void get_poc(l2str *l2rec, l2prodstr *p, float prod[]);
float poc_stramski_hybrid(float *Rrs, int32_t sensorID);
void get_fsat(l2str *l2rec, float flh[]);
void get_fsat2(l2str *l2rec, float flh[]);
void get_fqy(l2str *l2rec, float fqy[]);
void get_fqy2(l2str *l2rec, float fqy[]);
void get_ipar(l2str *l2rec, float ipar[]);
void get_par_scalar(l2str *l2rec, float par[]);
void get_par_below_surface(l2str *l2rec, float par[]);
void get_mu_cosine(l2str *l2rec, float mu[]);
void get_ipar2(l2str *l2rec, float ipar[]);
void get_ipar_below_surface(l2str *l2rec, float ipar[]);
void get_depth_classification(l2str *l2rec, float depth[]);
void get_par(l2str *l2rec, float par[]);
void get_par2(l2str *l2rec, float par[]);
void get_ipar_scalar(l2str *l2rec, float ipar[]);
void get_taucld (l2str *l2rec, float taucld[]);
void get_clfr (l2str *l2rec, float clfr[]);
void get_bsi(l2str *l2rec, float BSi[]);
void get_ssn(l2str *l2rec, float ssn[]);
void get_angstrom(l2str *l2rec, int band, float angst[]);
void get_ms_epsilon(l2str *l2rec, float eps[]);
void get_es(l2str *l2rec, int band, float es[]);
void get_toa_refl(l2str *l2rec, int band, float rhot[]);
void get_dust_index(l2str *l2rec, float dust[]);
void get_ndvi(l1str *l1rec, float ndvi[]);
void get_evi(l1str *l1rec, float evi[]);
void get_evi2(l1str *l1rec, float evi2[]);
void get_evi3(l1str *l1rec, float evi3[]);
void get_cci(l1str *l1rec, float cci[]);
void get_ndii(l1str *l1rec, float ndii[]);
void get_ndsi(l1str *l1rec, float ndsi[]);
void get_ndwi(l1str *l1rec, float ndwi[]);
int get_hyper_vi(l1str *l1rec, int product_number, float product[]);
void get_smoke(l2str *l2rec, float smoke[]);
void get_Kd(l2str *l2rec, l2prodstr *p, float Kd[]);
void get_photic_depth(l2str *l2rec, l2prodstr *p, float Z[]);
void cdom_morel(l2str *l2rec, l2prodstr *p, float prod[]);
void cdom_mannino(l2str *l2rec, int prodnum, float prod[]);
  //void get_soa(l2str *l2rec, int32_t prodID, float prod[]);
  //int run_soa_sma(l2str *l2rec, int32_t ip);
void vcal(l2str *l2rec, l2prodstr *p, float vcal[]);
float aw_spectra(int32_t wl, int32_t width);
float bbw_spectra(int32_t wl, int32_t width);
void get_aw_bbw(l2str *l2rec, float wave[], int nwave, float *aw, float *bbw);
float seawater_nsw(float wave, float sst, float sss, float *dnswds);
float seawater_bb(float wave, float sst, float sss, double delta);
void seawater_set(l1str *l1rec);
float seawater_get_n(int32_t ip, int32_t ib);
float seawater_get_a(int32_t ip, int32_t ib);
float seawater_get_bb(int32_t ip, int32_t ib);
void get_bbws(l2str *l2rec, l2prodstr *p, float prod[]);
void get_avw(l2str *l2rec, float avw[]);
void get_Rrs_brightness(l2str *l2rec, float Rrs_brightness[]);
void get_lambda_max(l2str *l2rec, float lambda_max[]);
void get_Cphyt(l2str *l2rec, float cphyt[]);
void get_Cpicophyt(l2str *l2rec, l2prodstr *p, float *cell_abundance);

int atmocor1_land(l1str *l1rec, int32_t ip);
void polcor(l1str *l1rec, int32_t ip);
int get_rhos(l1str *l1rec, int32_t ip);
int8_t *get_qual_sst(l2str *l2rec);
int8_t *get_qual_sst4(l2str *l2rec);
int8_t *get_qual_sst_triple(l2str *l2rec);
int16_t *get_flags_sst(l2str *l2rec);
int16_t *get_flags_sst4(l2str *l2rec);
int16_t *get_flags_sst_triple(l2str *l2rec);
float *get_sst_dust_correction(l2str *l2rec);
float *get_sst(l2str *l2rec);
float *get_sst4(l2str *l2rec);
float *get_sst_triple(l2str *l2rec);
float *get_bias_sst(l2str *l2rec);
float *get_bias_sst4(l2str *l2rec);
float *get_bias_sst_triple(l2str *l2rec);
float *get_stdv_sst(l2str *l2rec);
float *get_stdv_sst4(l2str *l2rec);
float *get_stdv_sst_triple(l2str *l2rec);
float *get_bias_mean_sst(l2str *l2rec);
float *get_bias_mean_sst4(l2str *l2rec);
float *get_bias_mean_sst_triple(l2str *l2rec);
int16_t *get_counts_sst(l2str *l2rec);
int16_t *get_counts_sst4(l2str *l2rec);
int16_t *get_counts_sst_triple(l2str *l2rec);

float get_sstref(short reftyp, char *file, l1str *l1rec, int32_t ip);
float get_sssref(char *file, float lon, float lat, int day);
void calcite(l2str *l2rec, l2prodstr *p, float caco3[]);
void tindx_morel(l2str *l2rec, int32_t ip, float *tindx);
void tindx_shi(l2str *l2rec, int32_t ip, float *tindx);
float conv_rrs_to_555(float Rrs, float wave, float dRrs_in, float *dRrs_out);

float water_vapor(int ib, float wv, float airmass);
int ice_mask_init(char *file, int year, int day, float threshold);
char ice_mask(float lon, float lat);
float ice_fraction(float lon, float lat);
void get_ice_frac(l2str *l2rec, float ice[]);
void get_tauc(l2str *l2rec, float tauc[]);
void get_mgiop(l2str *l2rec, l2prodstr *p, float prod[]);
void get_gsm(l2str *l2rec, l2prodstr *p, float prod[]);
int16_t *get_iter_gsm(l2str *l2rec);
void iops_gsm(l2str *l2rec);
void get_giop(l2str *l2rec, l2prodstr *p, float prod[]);
int16_t *get_iter_giop(l2str *l2rec);
int16_t *get_flags_giop(l2str *l2rec);
void iops_giop(l2str *l2rec);
void get_carder(l2str *l2rec, l2prodstr *p, float prod[]);
int16_t *get_flags_carder(l2str *l2rec);
void iops_carder(l2str *l2rec);
void chl_carder_empirical(l2str *l2rec, float prod[]);
void get_pml(l2str *l2rec, l2prodstr *p, float prod[]);
void iops_pml(l2str *l2rec);
void get_qaa(l2str *l2rec, l2prodstr *p, float prod[]);
unsigned char *get_flags_qaa(l2str *l2rec);
void iops_qaa(l2str *l2rec);
void get_niwa(l2str *l2rec, l2prodstr *p, float prod[]);
void iops_niwa(l2str *l2rec);
int16_t *get_flags_niwa(l2str *l2rec);
void iops_las(l2str *l2rec);
int get_bbp_qaa(l2str *l2rec, int ip, float tab_wave[], float tab_bbp[], int tab_nwave);
int get_bbp_las(l2str *l2rec, int ip, float tab_wave[], float tab_bbp[], int tab_nwave);
float get_bbp_las_eta(l2str *l2rec, int ip);
float get_bbp_las_eta_ksm(l2str *l2rec, int ip);
void get_iops(l2str *l2rec, int32_t iop_opt);
void set_iop_flag(float wave[], int32_t nwave,
        float a[], float aph[], float adg[],
        float bb[], float bbp[], int16_t *flag);
float aph_bricaud(float wave, float chl);
float aph_ciotti(float wave, float sf);
float get_aphstar(float wave, int dwave, int ftype, float proxy);
float rrs_above_to_below(float Rrs);
void optical_class(l2str *l2rec, l2prodstr *p, float prod[]);
float *get_class_ward_owmc(l2str *l2rec);
float *get_class_k_owmc(l2str *l2rec);
float *get_class_34k_w_owmc(l2str *l2rec);

  /*
void myprod1(l2str *l2rec, float prod[]);
void myprod2(l2str *l2rec, float prod[]);
void myprod3(l2str *l2rec, float prod[]);
void myprod4(l2str *l2rec, float prod[]);
void myprod5(l2str *l2rec, float prod[]);
void myprod6(l2str *l2rec, float prod[]);
void myprod7(l2str *l2rec, float prod[]);
void myprod8(l2str *l2rec, float prod[]);
void myprod9(l2str *l2rec, float prod[]);
void myprod10(l2str *l2rec, float prod[]);
  */

float westernmost(float lon1, float lon2);
float easternmost(float lon1, float lon2);

/* Filter functions */

void fctl_init(fctlstr *fctl);
int fctl_set(fctlstr *fctl, int32_t npix, char *fname,
        int32_t band, int32_t nx, int32_t ny, int32_t minfill, int32_t nbands);
void filter(fctlstr *fctl, l1qstr *l1que, l1str *l1rec, int32_t dscan);
int rdfilter(char *file, fctlstr *fctl, int32_t nbands);
void fdilate(l1qstr *l1que, int32_t nx, int32_t ny, int flag, char kernel[], l1str *l1rec);
void fclean(l1qstr *l1que, int32_t nx, int32_t ny, int flag, char kernel[], l1str *l1rec);
void fLTmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fLTRmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fLTmed(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fLTRmed(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fEPSmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fEPSiqmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fLTRiqmean(l1qstr *l1que, int32_t nx, int32_t ny, int ib, int32_t minfill, char kernel[], l1str *l1rec);
void fLTrat(l1qstr *l1que, int32_t nx, int32_t ny, l1str *l1rec);

/* viirs pixel conversion */
void viirs_pxcvt_2uag(int, int *, int *);
void viirs_pxcvt_2ag(int, int *);
void viirs_pxcvt_agdel(int, int, int *);

/* Fortran functions called from C */

void sunangs_(int *iyr, int *iday, float *gmt, float *xlon, float *ylat,
        float *sunz, float *suna);
int getglint_(float *, float *, float *, float *, float *, float *);
int getglint_iqu_(float *, float *, float *, float *, float *, float *, float *, float *);
void get_um_prod_(int32_t* idProd, int32_t* pix, float* val);
void clear_um_prod_();
  /*
int atmcor_soa_(
        int32_t *sensorID,
        char *sensorNm,
        int32_t *nwave,
        float *wave,
        int32_t *scan,
        int32_t *pixel,
        float *solz,
        float *senz,
        float *raz,
        float *lat,
        float *lon,
        float *Lt,
        float *rho_r,
        float *Fo,
        float *Tau_r,
        float *aw,
        float *bbw,
        float *aphstar,
        float *adg_s,
        float *bbp_s,
        float *Rs,
        float *Rw,
        float *wv,
        float *t_sen,
        float *t_sol,
        float *optTaua,
        float *optW0,
        float *optChl,
        float *optAcdm,
        float *pcentCDM,
        float *optBbp,
        float *optMr,
        float *optMi,
        float *optV,
        int32_t *status);
  */
  
void get_fdiff(l2str *l2rec, float fdiff[]);

void get_cdom_morel(l2str *l2rec, l2prodstr *p, float prod[]);

void optical_water_type(l2str *l2rec, l2prodstr *p, void *vptr);

void* prodgen(l2prodstr *p, l2str *l2rec);

void virtual_constellation(l2str *l2rec, l2prodstr *p, float prod[]);
void bioOptBandShift(l2str *l2rec, l2prodstr *p, float prod[]);

void get_swim(l2str *l2rec, l2prodstr *p, float prod[]);
void iops_swim(l2str *l2rec);

int compfloat(float *x, float *y);
void elev_init(char* elevFilename, char* elevAuxFilename);
float get_elev(float lat, float lon);

int read_target_l3(filehandle *file, l1str *l1rec, int32_t nbands, tgstr *tgrec);
int ncio_dim_siz(int, char *);
int ncio_grab_f_ds(int, char *, float *);
int ncio_grab_stdsclf_ds(int, char *, float, float *);

float get_mld(char* mldfile, float lon, float lat, int day);
void get_pft_hirata(l2str *l2rec, l2prodstr *p, float prod[]);
void get_pft_uitz(l2str *l2rec, l2prodstr *p, float prod[]);
void get_npp(l2str *l2rec, int32_t, float prod[]);
float chl_abi(l2str *l2rec, float nLw[]);

void run_raman_cor(l2str *l2rec, int ip);
void get_bpar(l2str *l2rec,  l2prodstr *p, float prod[]) ;

void get_habs_ci(l2str *l2rec, l2prodstr *p, float ci[]);
void get_habs_mph(l2str *l2rec, l2prodstr *p, float mph_chl[]);
uint8_t* get_flags_habs_mph(l2str *l2rec);
uint8_t* get_flags_habs(l2str *l2rec);
void get_psd_ksm(l2str *l2rec, l2prodstr *p, float prod[]);

int get_cmp( l2str *l2rec, int prodnum, float prod[]);
unsigned char *get_cmp_byt( l2str *l2rec, int32_t );

float* giop_get_aph_pointer();
float* giop_get_adg_pointer();
float* giop_get_bbp_pointer();
float* giop_get_bbp_s_pointer();
float* giop_get_ubbp_s_pointer();
float* giop_get_a_unc_pointer();
float* giop_get_bb_unc_pointer();
void run_giop(l2str *l2rec);
int giop_ran(int recnum);

float first_deriv(float x[], float y[], int n);

void gas_trans_uncertainty(l1str *l1rec);
int noise_index_oci(float wave);
void noise_model_ocis(l1str *l1rec, float *noise);
void noise_model_hmodisa(l1str *l1rec, float *noise);
void noise_model_hmodist(l1str *l1rec, float *noise);
void lt_agregat_ocis(l1str*l1rec);
float *get_uncertainty(l1str *l1rec);

void get_unc_fsat(l2str *l2rec, float uflh[]);
float get_aphstar_pderiv(float wave, int dwave, int ftype, float proxy);
void get_Cphyt_unc(l2str *l2rec, float cphyt_unc[]);

void get_sdp(l2str *l2rec, l2prodstr *p, float prod[]);

#ifdef __cplusplus
}
#endif

#endif
