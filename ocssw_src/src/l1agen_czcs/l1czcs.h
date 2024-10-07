#ifndef __CZCS_DEFINED
#define __CZCS_DEFINED 1
/*
Name: l1czcs.h

This file contains structures and typedef's for czcs lever 1 data.
Also, function prototypes are defined 

The typedefs:

     HEADER2_TYPE  :  typedef for a czcs level 1 tape
                      document record;
     DATA_REC_TYPE :  typdef for a czcs level 1 data record.

JPB 20 Jun 95

W. Robinson, SAIC, 5 Mar 2004   adapt for use in l1czcs program, to 
                                update for Watson Gregg's Algorithm
                                Improvement

 */

#define NCZCS_PIX 1968
/*
 *   definition of document record:
 */
#include "hdf.h"
#include "mfhdf.h"

struct bit_field_1 /* a bit field used in document record */ {
#ifdef __DEC_MACHINE        
    unsigned p_rec_num : 12;
    unsigned sp1 : 4;
    unsigned rec_id : 6;
    unsigned file : 2;
    unsigned valid_data_flg : 8;
#else
    unsigned p_rec_num : 12;
    unsigned sp1 : 4;
    unsigned file : 2;
    unsigned rec_id : 6;
    unsigned valid_data_flg : 8;
#endif
};

typedef struct bit_field_1 BIT_FLD_1_TYPE;

struct line_long {
    int slope;
    int intercept;
};

typedef struct line_long LINE_TYPE_LONG;

struct line_short {
    short slope;
    short intercept;
};

typedef struct line_short LINE_TYPE_SHORT;

struct hdr2 {
    BIT_FLD_1_TYPE file_info;

    unsigned char target_area_code[3];
    unsigned char file_no;

    int tape_seq_no;
    int film_frame_no;
    short starting_yr;

    short starting_day;

    int start_msec;
    int inc_in_msec;
    short orbit;
    short no_of_scans;
    short lat_cntr;
    short long_cntr;
    short lat_crnr_fitl;
    short long_crnr_fitl;
    short lat_crnr_fitr;
    short long_crnr_fitr;
    short lat_crnr_litl;
    short long_crnr_litl;
    short lat_crnr_litr;
    short long_crnr_litr;
    unsigned char ilt_flags;
    unsigned char parameter_presence;
    short no_of_missing_scans_all;
    short no_of_missing_scans[6];
    unsigned char alg_id_chan[6];
    unsigned char alg_id_location;
    char spare2;
    int decom_run_no;
    int decom_reel_no;
    short no_of_hdt_sync_losses;
    short no_of_hdt_parity_errs;
    short no_of_wbvt_sync_losses;
    short no_of_wbvt_bit_slips;
    short ave_subcomputed_data[32];
    char spare3;
    unsigned char bp_flag;
    short baseplate_temp;
    int spare4[134];
    unsigned char gain;
    unsigned char threshold;
    short tilt;
    short scene_cntr_yr;
    short scene_cntr_doy;
    int scene_cntr_msec;
    short solar_el_cntr;
    short solar_az_cntr;
    short cntr_roll;
    short cntr_pitch;
    short cntr_yaw;
    unsigned char top_bot_tick_label_flag;
    unsigned char lft_rgt_tick_label_flag;
    short top_lft_tick_label;
    short top_rgt_tick_label;
    short bot_lft_tick_label;
    short bot_rgt_tick_label;
    short lft_top_tick_label;
    short lft_bot_tick_label;
    short rgt_top_tick_label;
    short rgt_bot_tick_label;
    unsigned char top_tick_inc;
    unsigned char bot_tick_inc;
    unsigned char lft_tick_inc;
    unsigned char rgt_tick_inc;
    short top_tick_loc[27];
    short bot_tick_loc[27];
    short lft_tick_loc[27];
    short rgt_tick_loc[27];
    LINE_TYPE_LONG chan_line[6];
    short conversion_tab_chan6[256];
    LINE_TYPE_SHORT enhancement_eqs[6];
    int spare5[2];
    int czcs_1lt[945];
    unsigned char dummy[304]; /* rounded out record to 512 blk boundary */
};

typedef struct hdr2 HEADER2_TYPE;

/*
 *      definition of data record:
 */



struct bit_field_2 { /* bit field used in data record */
#ifdef __DEC_MACHINE
    unsigned p_rec_num : 12;
    unsigned sp1 : 4;
    unsigned rec_id : 6;
    unsigned file : 2;
    unsigned undefined : 3;
    unsigned volt_outside_range : 1;
    unsigned cal_outside_range : 1;
    unsigned channel_not_present : 1;
    unsigned qstnble_attitude : 1;
    unsigned qstnble_ephemeris : 1;
#else
    unsigned p_rec_num : 12;
    unsigned sp1 : 4;
    unsigned file : 2;
    unsigned rec_id : 6;
    unsigned qstnble_ephemeris : 1;
    unsigned qstnble_attitude : 1;
    unsigned channel_not_present : 1;
    unsigned cal_outside_range : 1;
    unsigned volt_outside_range : 1;
    unsigned undefined : 3;
#endif
};

typedef struct bit_field_2 BIT_FLD_2_TYPE;

struct voltage_step { /*  a component of voltage staircase */
    unsigned char whole;
    unsigned char fraction;
};

typedef struct voltage_step VOLTAGE_STEP;

struct data_record_struct {
    BIT_FLD_2_TYPE info;
    short seq_no;
    unsigned char spare1;
    unsigned char time_update_flag;
    short year;
    short doy;
    int msec_of_day;
    short subcommuted_data_val_cnt;
    unsigned char subcom_id;
    unsigned char spare2;
    VOLTAGE_STEP volt_stair[6][16]; /* one 16 step staircase per channel */
    short cal_lamp_rad[6];
    short blk_bdy_tmp_cnt;
    short wbvt_bit_slips_summary;
    short hdt_sync_losses;
    short hdt_parity_errs;
    short wbvt_bit_slips;
    int lat_anchr_pts[77];
    int long_anchr_pts[77];
    unsigned short pxl_no_nadir;
    unsigned char cal_qual[6]; /* one per channel */
    unsigned char radiance_cnts[6][1968]; /* 1968 counts per channel */
    int spare3[28];
    unsigned char dummy[20]; /* round out record to 512 block boundary */
};

typedef struct data_record_struct DATA_REC_TYPE;

/*
 * Definition of VAX header rec :
 */

struct VAX_header_struct {
    unsigned short magic1;
    unsigned short magic2;
    unsigned short data_rec_len;
    unsigned short num_of_doc_recs;
    unsigned short offset_1st_data_rec;
    unsigned short type_code;
    unsigned short num_of_data_recs;
    unsigned short orbit;
    unsigned short year;
    unsigned short header_rec_offset;
    unsigned short header_rec_len;
    unsigned short doc_rec_len;
    unsigned short dummy1;
    unsigned short dummy2;
    unsigned short dummy3;
    unsigned short scanner_tilt;
    char text[480];
};

typedef struct VAX_header_struct HEADER1_TYPE;

/*
 *  structure to contain l1 and lat, lon anchor point data
 */

struct l1_data_struc_def {
    unsigned char *counts[6]; /* counts for 6 bands in form pix x line x band */
    /* with pixels running fastest  */
    int *msec;
    float *ctl_pt_lat;
    float *ctl_pt_lon;
    int *ctl_pt_cols; /* control pint columns and rows, evenly spaced   */
    int *ctl_pt_rows; /* and lastpoint on the line   */
    float *tilt; /* tilt / line */
    float *slat; /* start, center, end lat, lon per line */
    float *slon;
    float *clat;
    float *clon;
    float *elat;
    float *elon;
    unsigned char *cal_sum; /* calibration quality summary and */
    unsigned char *cal_scan; /* quality per scan */
    float *orb_vec;
    float *att_ang;
    float *pos_err;
    float *slope;
    float *intercept;
    short *gain;
#ifdef GEOM_CAL
    float *sen_zen; /* geometry and calibrated radiances for test */
    float *sen_az;
    float *sol_zen;
    float *sol_az;
    float *all_lat;
    float *all_lon;
    float *Lt_443;
    float *Lt_520;
    float *Lt_550;
    float *Lt_670;
    float *Lt_750;
    float *Lt_11500;
#endif
};

typedef struct l1_data_struc_def l1_data_struc;

/*
 *  Global attributes for output to the hdf file
 */

struct gattr_struc_def {
    char prod_name[40]; /* Product Name */
    char f_title[40]; /* Title */
    char datatype[16]; /* Data Type */
    char datacenter[256]; /* Data Center */
    char stn_name[256]; /* Station Name */
    float stn_lat; /* Station Latitude */
    float stn_lon; /* Station Longitude */
    char mission[40]; /* Mission */
    char mission_char[256]; /* Mission Characteristics */
    char sensor[40]; /* Sensor */
    char sensor_char[256]; /* Sensor Characteristics */
    char repl_flg[40]; /* Replacement Flag */
    char soft_id[128]; /* Software ID */
    char process_time[17]; /* Processing Time */
    char input_files[1280]; /* Input Files */
    char proc_ctl[1024]; /* Processing Controls */
    char start_time[17]; /* Start Time */
    char end_time[17]; /* End Time */
    char center_time[17]; /* Scene Center Time */
    short start_year; /* Start Year */
    short start_day; /* Start Day */
    int start_msec; /* Start Milisec */
    short end_year; /* End Year */
    short end_day; /* End Day */
    int end_msec; /* End Milisec */
    char start_node[11]; /* Start Node */
    char end_node[11]; /* End Node */
    int orbit; /* Orbit Number */
    int gain; /* Scene Gain */
    int thresh; /* Thresh */
    float tilt; /* Sensor Tilt */
    char lat_unit[15]; /* Latitude Units */
    char lon_unit[15]; /* Longitude Units */
    float center_lat; /* Scene Center Latitude */
    float center_lon; /* Scene Center Longitude */
    float cntr_sol_zen; /* Scene Center Solar Zenith */
    float up_lft_lat; /* Upper Left Latitude */
    float up_lft_lon; /* Upper Left Longitude */
    float lo_lft_lat; /* Lower Left Latitude */
    float lo_lft_lon; /* Lower Left Longitude */
    float up_rgt_lat; /* Upper Right Latitude */
    float up_rgt_lon; /* Upper Right Longitude */
    float lo_rgt_lat; /* Lower Right Latitude */
    float lo_rgt_lon; /* Lower Right Longitude */
    float start_cntr_lat; /* Start Center Latitude */
    float start_cntr_lon; /* Start Center Longitude */
    float end_cntr_lat; /* End Center Latitude */
    float end_cntr_lon; /* End Center Longitude */
    int pix_per_scan; /* Pixels per Scan Line */
    int scan_lines; /* Number of Scan Lines */
    int lac_pixl_start_no; /* LAC Pixel Start Number */
    int lac_pixl_subsample; /* LAC Pixel Subsampling */
    int cntr_scn_line; /* Scene Center Scan Line */
    int filled_lines; /* Filled Scan Lines */
    float limits[4]; /* Northernmost Latitude, Southernmost 
                    Latitude, Westernmost Longitude, Easternmost Longitude */
    float slope[6]; /* calibration slope and intercept from 2nd */
    float intercept[6]; /* header = derived from voltage staircase */
    /* Calibration Slope and Calibration Intercept */
    float roll; /* Center Roll S/C roll, pitch and yaw at scene center */
    float pitch; /* Center Pitch */
    float yaw; /* Center Yaw */
    int ctl_pt_incr;
    int n_ctl_pt; /* Number of Pixel Control Points */
    int n_ctl_lin; /* Number of Scan Control Points */
    unsigned char ilt_flags; /* ILT Flags */
    unsigned char parm_presence; /* Parameter Presence Code */
    short n_miss_scans; /* Number of Missing Scan Lines */
    short n_scan_mis_chan[6]; /* Number of Scans with Missing Channels */
    short n_hdt_sync_loss; /* Number of HDT Sync Losses */
    short n_hdt_parity_err; /* Number of HDT Parity Errors */
    short n_wbvt_sync_loss; /* Number of WBVT Sync Losses */
    short n_wbvt_slips; /* Number of WBVT Slip Occurrences */
};

typedef struct gattr_struc_def gattr_struc;

struct timqual_struc_d {
    int nscan; /* # scan lines */
    int n_ctl_pt; /* # pixel ctl points - needed to set up out file */
    int *qual; /* quality - 0 good, 1 bad */
    int32_t *msec; /* msec from the file */
};

typedef struct timqual_struc_d timqual_struc;

struct mstr_struc_d {
    int32_t *msec; /* time tag */
    short *exist; /* does this scan exist? 1 = yes */
    short *qual; /* quality 0 = good, 1 bad */
    short *ds_num; /* dataset # in file list to take scan from */
    int32_t *out_scan; /* output scan line to place the data */
    short *scan; /* scan # to use */
    int32_t *in_msec; /* the in_ values are for new candidate scans */
    short *in_exist;
    short *in_qual;
    short *in_scan;
};

typedef struct mstr_struc_d mstr_struc;

#endif

/*
 *  the function prototypes
 */
int main(int, char *[]);
int read_crtt(char *, gattr_struc *, l1_data_struc *);
void get_record_info(HEADER1_TYPE, short *, short *, short *, short *,
        short *);
short reverse_short_int(short);
int32_t reverse_long_int(int32_t);
void hdr_2_gattr(HEADER2_TYPE, gattr_struc *);
float fixed_pt_2_floating_pt(int, int);
int czcs_l1_write(char *, l1_data_struc, gattr_struc);
int create_global_attribute(char *, int, gattr_struc);
int create_band_sds(int, int, unsigned char *[], int, int);
void czcs_ctl_pt(DATA_REC_TYPE, gattr_struc *, int, l1_data_struc *);
void lonlat_(float[], float[], float[], float[]);
void lladjust_(float *, float *, float *, float *, float[], float[], float[]);
void satang_(double *, double *, float *, float *, float *, float *, float *,
        float *, float *, float *);
void sunangs_(int *, int *, float *, float *, float *, float *, float *);
int time_str(short, short, int, char *);
void cz_ll_upd(l1_data_struc *, gattr_struc *);
int wrt_czcs_sla(int32, int32, int, l1_data_struc);
int wrt_czcs_qual(int32, int32, int, l1_data_struc);
int32 set_czcs_ctl_data(int32, int32, gattr_struc, l1_data_struc);
int cz_clean(gattr_struc *, l1_data_struc *);
/* new for the merge program */
int read_file_list(char *, char **, int);
void usage(char *);
void olap_resolve(mstr_struc *, int, int, int, int);
void fill_mstr(int *, mstr_struc *, timqual_struc *, int, int, int);
int cztimqual(char *, timqual_struc *, int *);
int cz_l1_read(char *, int, gattr_struc *, l1_data_struc *);
int cz_dat_alloc(int, int, int, l1_data_struc *);
void cz_dat_free(l1_data_struc *, int);
int cz_mov_scn(int, int, char *, mstr_struc *, int, int, gattr_struc *,
        l1_data_struc *);
void cz_meta_adj(l1_data_struc *, gattr_struc *);
void cz_sd_set(l1_data_struc *, gattr_struc *);
int fill_orb_dat(l1_data_struc *l1_data, gattr_struc *gattr);
#ifdef __DEC_MACHINE
void convert_data_rec_to_dec(DATA_REC_TYPE *);
void convert_header_rec_to_dec(HEADER2_TYPE *);
#endif
