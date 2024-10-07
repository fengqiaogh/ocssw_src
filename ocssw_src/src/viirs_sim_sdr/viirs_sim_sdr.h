/*******************************************************************

   viirs_sim_sdr.h

   purpose: include file for the viirs sdr simulation program

   Parameters: 
      Type              Name            I/O     Description
      ----              ----            ---     -----------

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson, SAIC 23-Sep-2008     Original development
      W. Robinson, SAIC 20 Nov 2008     modify for line-by-line I/O,
                                        structure seperation for in, out, 
                                        and info

 *******************************************************************/

/*
 *  Note that hdf5.h is needed for this
 */
#include "h5io.h"
#include "readL2scan.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

#define NDET 16
#define N_VNIR_BND 7
#define N_IR_BND 9
/* assume 1 geo file */
#define MAX_BND N_VNIR_BND + N_IR_BND
#define MAX_FILES MAX_BND + 1
#define RAD_CGS_2_MKS 10.
#define SEC_PER_SCAN 1.7864
#define PI 3.141592654
/*
 *  some needed missing data fill value indicators
 */
#define ONBOARD_PT_FLOAT32_FILL  -999.7
#define ONBOARD_PT_UINT16_FILL  65533
#define SOUB_UINT16_FILL 65528
#define ERR_UINT16_FILL 65531
#define MISS_UINT16_FILL 65534
#define MISS_FLOAT32_FILL -999.8
/*
 *  VIIRS scan angle types
 */
#define VIR_SCAN_AGSMP 0   /* aggregated sample # */
#define VIR_SCAN_UASMP 1   /* unaggregated sample # */
#define VIR_SCAN_ANG 2     /* scan angle, degrees */
#define VIR_SCAN_AOI 3     /* AOI on HAM mirror, degrees */

/*
 *  input controls to program
 */
struct ctl_struc_d {
    char in_geo_file[500]; /*  1st arg, the geo data file with basic 
                              data in it */
    char out_loc[500]; /*  output SDR location.  Files will be named
                          VYYYYDDDHHMMSS_GEO.h5 for geo and
                          VYYYYDDDHHMMSS_SDRMXX.h5 for the band files,
                          XX is the band number  */
    int l2_use; /* flag for existance (1) or absence (0) of an L2 file with
                  TOA radiances needed for the SDR */
    char l2_file[500]; /* name of the input L2 file */
    int meta_use; /* 1 if metadata output is desired, otherwise 0 */
    char meta_file[500]; /* name of file to put metadata in */
    int rhos_use; /* flag for existance (1) or absence (0) of a surface 
                    reflectance file for filling gaps in TOA data*/
    char rhos_file[500]; /* name of surface reflectance file */
    int rhos_opt; /* reflectance data use option: 0 - replace where Lt 
                     from ocean color is missing, 1 - replace everywhere */
    int bowtie_opt; /* include the bow tie deletion if 1, not if 0 */
    int any_artifact; /* are there any artifacts? 0 = no, 1 = yes */
    int out_scn_fmt; /* switch to modify scan format: 2 = default - do not 
      change input data format, 1 - if unaggregated with margin, remove the 
      margin, 0 - make output aggregated = standard SDR format */
    int oxt_mode; /* optical crosstalk (OXT) mode: 0 (default) no Xtalk applied,
                    1 apply it */
    char oxt_coef[500]; /* OXT influence coefficients, default is 
                    $VIIRS_SIM_DATA/oxt_default.h5 */
    char inter_band[500]; /* Inter-band data table, default is 
                       $VIIRS_SIM_DATA/inter_band_default.dat */
    int make_m; /* control to make M bands: 0 - do the NIV NIR 1-7, 
                1 - do all 16 bands */
    int sdr_overwrite; /* control to allow overwrite of an existing SDR file:
                0 - (default) do not allow, 1 - allow overwrite */
    /*  for removing vicarious gain and offset from the TOA from l2gen_inv 
        if ever needed */
    int vic_cal_chg; /* 0 if no change, 1 if a change in gain or offset */
    /* from default of gain 1, offset 0 */
    float gain[MAX_BND]; /* gains  */
    float offset[MAX_BND]; /* offset, eq: Lt = ( Lt0 * gain ) + offset  */
    /* create time control */
    char cre_time[23]; /* create time in form YYYYMMDD/HHMMSS.FFFFFF */
    int fname_opt; /* file name option: 0 to make std name, 
                     1 to omit create time */
    /*  for the count conversion and electronic crosstalk  */
    int count_cal_opt; /* cal option: 0 - no count calibration, 1 - do cal
                         and de-cal (set if electronic Xtalk is set),
                         2 - as (1), but also integerize the dn value */
    char count_cal_gain_file[500]; /* file containing the gain and dark
                                     for count to radiance conversion */
    char count_cal_rvs_file[500]; /* file containing the RVS for count
                                    to radiance conversion */
    char count_decal_gain_file[500]; /* as above, but for de-calibration
                                       if different from the cal */
    char count_decal_rvs_file[500]; /* as count_cal_rvs_file, but for de-cal */
    int count_dark_opt; /* set to 0 for no dark subtraction, 1 to do subtract */
    int ext_opt; /* electronic crosstalk mode: 0(default) - no crosstalk, */
    /* 1 - perform the crosstalk  */
    char ext_coeff_file[500]; /* name of ext file or 'Unspecified' for */
    /* use of trivial coefficients (zero values) */
    char id_domain[500]; /* domain designator in dataset and attributes */
    char id_origin[500]; /* origin designator in dataset and attributes */
    /*  noise controls - an option and a file of coeffs */
    int noise_mode; /* 0 if not adding noise, 1 if adding it */
    char noise_coef[500]; /* noise coeff file */
    /*  stray light controls - an option and a file of coeffs */
    int stray_opt; /* stray light use option 0 - no, 1 - yes */
    char stray_tbl[500]; /* stray light PSF coeff file */
};

typedef struct ctl_struc_d ctl_struc;

typedef long long int64;
typedef unsigned long long uint64;

/*
 *  general data information from the simulated geo file and synthesized L2
 */
struct sdr_info_struc_d {
    char origin[50]; /* Originator of the data, for file name create and attr */
    char domain[50]; /* Domain - for file name create */
    int year; /*  information about the geo data - 4 digit year*/
    int day; /*  day of year */
    double start_sec; /*  data start second */
    char ofile_base[200]; /* storage for the output file base with time */
    char sdr_files[MAX_FILES][200]; /* fil names of SDR files */
    /*  Some derived values that are convienient to carry here */
    int64 t58_day; /* time of start of the day in microsec past 1/1/1958 */
    char cre_date[10]; /* string of the create data of file: YYYYMMDD */
    char cre_time[20]; /* string of the create time of file: HHMMSS.SSSSSSZ */
    int64 st_58_t; /* data start in micro sec > 1958 */
    int64 en_58_t; /* data end in micro sec > 1958 */
    char st_date[10]; /* data start as variable cre... (cre_date... above) */
    char st_time[20]; /* data start time as cre... */
    char en_date[10]; /* data end as cre... */
    char en_time[20]; /* data end time as cre... */
    float *geo_pos; /* position data from geolocation file, in km */
    float *geo_vel; /* velocity data from geolocation file, in km/s */
    float *geo_att; /* sensor attitude data from geolocation file, 
                         in degrees */
    unsigned char *ham_side; /* HAM half angle mirror side storage,
                      just set to 0 on 1st scan 1 on next ... */
    double *scan_time; /* read from the geo file, scan time in sec of day */
};

typedef struct sdr_info_struc_d sdr_info_struc;

struct in_rec_struc_d {
    int year; /* year of data */
    int yday; /* day of year */
    int msec; /* millisec of day for the scan */
    int nbnd; /* # bands */
    int npix; /* # pixels  */
    int nlin; /* # lines in whole file */
    int nscan; /* # scans should be nlin / 16 */
    /*  input format info read from the geo file */
    int scn_fmt; /* input scan format, found from geo file: 0 aggregated,
                      1 unaggregated - 6304 pix, 2 unaggregated with margin */
    int margin[2]; /* along-track and -scan margins around std scan */
    int ndet_scan; /* convenience - NDET + 2 * margin[0] = # det in a scan */
    unsigned char ham_side; /* HAM half angle mirror side for the scan */
    /*  */
    h5io_str geo_fid; /*  file id for geo file */
    h5io_str geo_dat_id[6]; /* dataset ids for the input geo info:
                             lat, lon, sena, senz, sola, solz   */
    float *lat; /*  the line[s] storage arrays  */
    float *lon;
    float *sena;
    float *senz;
    float *sola;
    float *solz;
    float *bnd_lt[MAX_BND]; /* for the radiances (reflectances made on
                           output) in W m^-2 sr^-1  or BBT in K */
    unsigned char *bnd_q[MAX_BND]; /* flags for bad values in radiances, */
    /*  see bnd_q in out_rec */
    float *dn[MAX_BND]; /* if needed, storage for the count data in float */
    unsigned char *gain_bit[MAX_BND]; /* seperate gain state bit for counts */
    char *dn_sat[MAX_BND]; /* indicator that the dn resulting from the radiance
                            would have been saturated in VIIRS if 1 */
    float ll_lims[6]; /*  These are latitude, longitude limits for the 
                         input geo data: low lat, high, low (-)lon, high
                         low (+)lon, high */
    int lam_band[MAX_BND]; /* the band wavelength in nm */
    float f0[MAX_BND]; /* the season corrected F0 value for the bands */
    /* info for the input L2 file with TOA rads (and atmos info for rhos work) */
    l2_prod l2_str; /*  The information structure for reading a l2 file */
};
typedef struct in_rec_struc_d in_rec_struc;

struct out_rec_struc_d {
    int year; /* year of data */
    int yday; /* day of year */
    int msec; /* millisec of day for the scan */
    int nbnd; /* # bands */
    int npix; /* # pixels  */
    int nlin; /* # lines  */
    int nscan; /* # scans should be nlin / 16 */
    /*  input format info read from the geo file */
    int scn_fmt; /* input scan format, found from geo file: 0 aggregated,
                      1 unaggregated - 6304 pix, 2 unaggregated with margin */
    int margin[2]; /* along-track and -scan margins around std scan */
    int ndet_scan; /* convenience - NDET + 2 * margin[0] = # det in a scan */
    /*    */
    h5io_str sdr_fid[MAX_FILES]; /*  file id for sdr files making output SDR */
    h5io_str sdr_dat_gid[2][MAX_FILES]; /* 2 groups under which lies the 
                                            geo or band datasets */
    h5io_str geo_dat_id[6]; /* dataset ids for the output geo info:
                             lat, lon, sena, senz, sola, solz   */
    h5io_str bnd_dat_id[2][MAX_BND]; /* dataset ids for the output band
                                           info: radiance, refl / BBT values */
    h5io_str qual1_m_id[MAX_BND]; /* QF1_VIIRSMBANDSDR dataset id */
    float *lat; /*  the line[s] storage arrays  */
    float *lon; /*  NOTE that for now, these will point to the input arrays */
    float *sena;
    float *senz;
    float *sola;
    float *solz;
    float *bnd_lt[MAX_BND]; /* for the radiances (reflectances made on
                           output) in W m^-2 sr^-1 or BBT in K */
    unsigned char *bnd_q[MAX_BND]; /* flags for bad values in radiances, 0 if */
    /*  none and 1 if bow tie deletion on this pixel */
    /*  2 if bad value resultant from rads and refl fill */
    unsigned char *qual1_m[MAX_BND]; /* This is the QF1_VIIRSMBANDSDR dataset
          information for the current scan.  For an aggregated output 
          (scn_fmt = 0) it is std def of bits 2 and 3 having saturatio:
          0 no samples saturated, 1 some, 2 all saturated.  For unaggregated,
          (scn_fmt = 1, 2) bit 2 is gain state: 0 hi, 1 low and bit 3 is 
          saturation: 0 unsat, 1 saturated */
    float *f0; /* the season corrected F0 value for the bands */
    int *lam_band; /* the band wavelength in microns = nm / 10^3 */
    float scale[MAX_BND]; /* scale and offset for the VIIRS output radiance */
    float offset[MAX_BND];
    float refl_scale[MAX_BND]; /* as above but for reflectivity (reflective 
                                 bands) or BBT (emissive bands) */
    float refl_offset[MAX_BND];
    char out_bnd_typ[MAX_BND]; /* 0 if this band stored as unsigned short, 1 
                               for float  */
    char meas_typ[MAX_BND]; /* measurement type: 0 - rad + refl, 
                             1 - BBT and rad */
};
typedef struct out_rec_struc_d out_rec_struc;

/*
 *  This is an attribute description storage structure for outputting the 
 *  many attributes
 */
struct h5attr_struc_d {
    int express; /* a 0 to not output this entry as a attribute, 1 to do it */
    int type; /*  type of data - 0 most types, 1 - char type  */
    char *name; /*  attrib name */
    hid_t typ; /*  type of data, like H5T_NATIVE_INT, note for type = 1,
                  this is not used */
    int str_len; /*  string length - if 1 string, this can be 0, but
                    for an array, length of the strings (all same length)
                    This is not used for type = 0  */
    int ndim; /*  # dimensis of data array, 0 if a scalar */
    int *dim_siz; /*  lengths of the dimensions, not used if ndim = 0 */
    void *data; /*  the actual data to write */
};
typedef struct h5attr_struc_d h5attr_struc;

/*
 *  The next 2 structures contain the calibration information: the gain
 *  and the RVS information
 */
struct vir_gain_struc_d {
    int nham; /* # Half angle Mirror sides (2, I hope) */
    int ndet; /* # detectors (16 for VIIRS)  */
    int ngain; /* # gain states (2, I hope) */
    int nbnd; /* # bands - should be 11 for all reflective M bands,
                 but could be 7 for just vis and NIR */
    float *c0; /* zero order gain coefficient */
    float *c1; /* first order gain coefficient */
    float *c2; /* 2nd order gain coefficient */
    /* the lt is computed from counts, dn, with
       lt = c0 + c1 * dn + c2 * dn^2  */
    float *dark; /* dark count value, the dark-subtracted counts, dn, are
                   related to the raw counts, DN, by
                   dn = DN - dark */
    /* for c* and dark, the storage is f( HAM, det, gain, band )
       with band the fastest varying... */
    char gain_ct[16];
    /* gain_ct is the # gain states for the m band - it will
       determine if dn conversion will try to use the low gain
       or just stay with the high gain (more likely to
       saturate) */
};
typedef struct vir_gain_struc_d vir_gain_struc;

struct vir_rvs_struc_d {
    int nham; /* # Half angle Mirror sides (2, I hope) */
    int ndet; /* # detectors (16 for VIIRS)  */
    int nbnd; /* # bands - should be 11 for all reflective M bands,
                 but could be 7 for just vis and NIR */
    float aoi_range[2]; /* AOI range coefficients are applicable for  */
    float *a0; /* zero order rvs coefficient*/
    float *a1; /* 1st order rvs coefficient*/
    float *a2; /* 2nd order rvs coefficient*/
    /* array storage: f( HAM, det, band ) 
       with band the fastest varying... */
    /* RVS = a0 + a1 * AOI + a2 * AOI^2 */
};
typedef struct vir_rvs_struc_d vir_rvs_struc;

struct vir_straylt_struc_d {
    int nbands; /* # bands included in stray light table */
    int nrec_det; /* # lines (= detectors) in stray light table */
    /* note that a convolution array exists for each 
       band and line of a standard scan */
    int nsamp; /* # pixels in the convolution array */
    int ndet; /* # lines (= detectors) in convcolution array */
    int csamp; /* (0 origin) center sample of the convolution array */
    float *psf; /* the set of convolution arrays for each band and 
                  detector, as a 4 d array of 
                  ( nsamp, ndet, nbands, nrec_det ) with nsamp the 
                  fastest dimension */
};
typedef struct vir_straylt_struc_d vir_straylt_struc;

/*
 *  prototypes
 */
int check_usage(int, char *[], ctl_struc *);
int rd_sim_init(ctl_struc *, sdr_info_struc *, in_rec_struc *);
int rd_geo_init(ctl_struc *, sdr_info_struc *, in_rec_struc *);
int init_sdr(ctl_struc *, sdr_info_struc *, in_rec_struc *, out_rec_struc *);
int init_sdr_top(int, sdr_info_struc *, out_rec_struc *);
int init_geo_data(sdr_info_struc *, in_rec_struc *, out_rec_struc *);
int init_bnd_data(int, sdr_info_struc *, in_rec_struc *, out_rec_struc *);
int init_sdr_dpattr(int, h5io_str *, sdr_info_struc *);
int init_sdr_agg(int, h5io_str *, sdr_info_struc *);
int init_sdr_gran(int, h5io_str *, sdr_info_struc *, out_rec_struc *);
int rd_sdr_scan(int, ctl_struc *, sdr_info_struc *, in_rec_struc *);
int rd_geo_scan(int, sdr_info_struc *, in_rec_struc *);
int gen_const_rad_scn(sdr_info_struc *, int, in_rec_struc *);
int wr_sdr_scan(int, out_rec_struc *);
int wr_geo_scan(int, out_rec_struc *);
int wr_bnd_scan(int, out_rec_struc *);
int fin_sdr(ctl_struc *, in_rec_struc *, out_rec_struc *);
int wr_attr_seq(h5io_str *, int, h5attr_struc *);
int gen_sdr_fname(int, char *, sdr_info_struc *, int, char *);
int rd_rhos_scan(char *, int, int, int, float **);
int rhos_to_lt(int, float **, in_rec_struc *, int, sdr_info_struc *);
int scan_cvt(in_rec_struc *, out_rec_struc *);
int ang_avg(float, float, float, float, float *, float *);
int mod_artifact(ctl_struc *, in_rec_struc *);
int viirs_oxt(ctl_struc *, in_rec_struc *);
int viirs_oxt_ib_read(char *, float *, float *);
int viirs_oxt_comp(ctl_struc *, in_rec_struc *);
float bbt_2_rad(float, float);
int vir_xf_scan(float, int, int, float *);
int vset_cal_gain(char *, vir_gain_struc *);
int vset_cal_rvs(char *, vir_rvs_struc *);
int viirs_cal(ctl_struc *, in_rec_struc *);
int viirs_decal(ctl_struc *, in_rec_struc *);
int viirs_sim_input(int, char *[], ctl_struc *);
int32_t rdsensorinfo(int32_t, int32_t, const char *, void **);
int bnd_ix_2_sen_info(char *, void *);
int viirs_ext(ctl_struc *, in_rec_struc *);
int viirs_noise(ctl_struc *, in_rec_struc *, int);
int viirs_noise_coef_rd(char *, float *, float *, float *, float *);
int viirs_straylt(ctl_struc *, in_rec_struc *, int);
int viirs_straylt_rd(char *, vir_straylt_struc *);
/*  */
int jd_c(int, int, int);
int day2mday(int, int, int *, int *);
