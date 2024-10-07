#include "l1io.h"
#include <mfhdf.h>

/*******************************************************************

   fmt_check.h

   purpose: include file to set up some of the structures used
            in the format checker program

   Modification history:
      Programmer        Date            Description of change
      ----------        ----            ---------------------
      W. Robinson       17-Feb-1995     Original development
      L. Kumar		19-Dec-1995     Modified dims of some char variables
                                        of u_data
      W. Robinson       13-Jun-1996     add prototypes for routines
      W. Robinson       17-Sep-1996     add capability to do 66k bytes
                                        of data instead of 3000 in u_data
      W. Robinson       30-Sep-1996     Upgrade to include min / max
                                        check ranges for the sds values
      W. Robinson       31-Oct-1996     Upgrade to recognize float64
      W. Robinson       15-Mar-2001     make changes for linux usage
      W. Robinson, SAIC 29 Mar 2005     update to fully table driven model
                                        and call I/O for sections from here
      W. Robinson, SAIC 3 Feb 2010      add multipliers for float uncertainty 
                                        in checking of attributes to a constant

 *******************************************************************/

/*
 *  The FMT_8BLEN variable is defined only for the linux so fmt_check
 *  gets built with smaller arrays
 */
#ifdef FCK_8LEN
#define FMT_8BLEN FCK_8LEN
#else
#define FMT_8BLEN 8251
#endif

/*
 *  define some codes for use
 *
 *  attribute value read keywords
 */
enum att_rd_code {
    ATT_RD_NOREAD, ATT_RD_READ_NOCK, ATT_RD_READ_ONE_VAL,
    ATT_RD_READ_INCLUSIVE
};

/*
 *  u_data allows the placing of short amounts of attribute data
 *  into the correctly formatted storage
 */

typedef union u_data_d {
    char chr[FMT_8BLEN * 8];
    int32 i32[FMT_8BLEN * 2];
    float32 f32[FMT_8BLEN * 2];
    int16 i16[FMT_8BLEN * 4];
    char i8[FMT_8BLEN * 8];
    unsigned char ui8[FMT_8BLEN * 8];
    float64 f64[FMT_8BLEN];
} u_data;

/*
 *  With the size of the u_data now, and the need for a similar
 *  structure for the SDS value range, introduce a shorter version
 *  for the sds 2-value range.
 */
typedef union s_data_d {
    int32 i32[2];
    float32 f32[2];
    int16 i16[4];
    char i8[8];
    unsigned char ui8[8];
    float64 f64[20];
} s_data;

/*
 * attr_str has all the information about an attribute and its value
 * if it is fixed
 */

typedef struct attr_str_d {
    int32 dim_index; /* index into dimension definitions, -1 if */
    /* no connection made  */
    char obj_nm[100]; /* Name of the science dataset or 'gbl'
                             if refering to a global attribute  */
    char access_nm[200]; /* This is the name used to access the value */
    char int_nm[200]; /* This is the expected name of the variable */
    int32 type; /* This is the type of the value returned */
    int32 count; /* This is the expected count of the variable, 
                             if <0, do not check the count */
    int32 read; /* This is the value read switch,
                             ATT_RD_NOREAD - do not read value,
                             ATT_RD_READ_NOCK - read value but do not check,
                             ATT_RD_READ_ONE_VAL - read the value and check 
                                 it with what is read into data below 
                             ATT_RD_READ_INCLUSIVE - make sure the value(s) are
                                 within the range described below  */
    u_data data; /* this is the expected value or range of the attribute */
} attr_str;

/*
 * sds_info_str will describe the science dataset and point to 
 * the attribute descriptions for the included attributes
 */

typedef struct sds_info_str_d {
    char name[100]; /* Name of the science dataset */
    int32 rank; /* # dimensions in this dataset */
    int32 type; /* type of the sds data */
    int32 e_ranges[3]; /* expected range of dimension >0 - check
                              with this value, 0 - do not check, -(value)
                              check vs the dimension size in the location:
                              'value' in the dimension id structure array */
    int32 n_attr; /* # of attributes in this (I have a max 
                              of 50 now) */
    int32 attr_ind[50]; /* place for index to attribute descriptions */
    int32 byt_per_val; /* # bytes in a value of this type  */
    int32 flg_rng; /* flag for min / max range checking: 0 - don't
                              check,  1 - do the check */
    s_data sds_rng; /* minimum and maximum (inclusive)for sds 
                              value checking */
} sds_info_str;

/*
 *  dim_id_str is the structure to hold all the dimension description 
 *  information.
 */
typedef struct dim_id_d {
    char att_nm[200]; /* attrinute name */
    char att_short[100]; /* short name for attr without spaces */
    int dim_size; /* the actual dimension size for this case */
} dim_id_str;

/*
 *  ras_str is the structure containing the raster definition(s)
 */
typedef struct ras_str_d {
    char lbl_ras[200]; /* raster label, *NONE* or no label */
    char lbl_pal[200]; /* palette label */
    int npix_indx; /* index into dimension struct array of the raster # pixels */
    int nlin_indx; /* index into dimension struct array of the raster # lines */
} ras_str;

/*
 * vg_info_str  describes the Vgroup products
 * The 3 Vgroups for bin description are there as well as the 
 * variable product descriptions
 */
typedef struct vg_info_str_d {
    int32 n_fields; /*  # Vdata fields in this Vdata */
    int32 num_typs[7]; /*  list of # types of each field */
    int32 nrec; /*  # Vdata records or -1 to use the nbin
                                 value passed in and only check that
                                 there are >= nbin # records in the Vdata  */
    char name[VSNAMELENMAX]; /*  name of this Vdata */
    char class[VSNAMELENMAX]; /*  class of this Vdata */
    char field_list[200]; /*  names of all the fields comma seperated */
} vg_info_str;

/*
 * structure for the entire format
 */
typedef struct fmt_str_d {
    int n_attr; /* # attribute descriptions */
    attr_str *att; /* attribute definitions */
    int n_sds; /* # SDS descriptions */
    sds_info_str *sds_info; /* SDS definitions */
    int n_dim_defs; /* # of dimension definitions */
    dim_id_str *dim_id; /* dimension definitions */
    int n_raster; /* # raster definitions (either 0 or 1 if a 
                              raster exists */
    ras_str *ras; /* raster definitions */
    int n_vgroup; /* # v group definitions (either 0 or 1 if a 
                              l3 set of v groups exist */
    int resolve; /* the resolve value about the bin size in km */
    int vg_nbin_indx; /* index of the '# of bins' in the dimension 
                  struct array designating the # of bins to be cross checked. */
    vg_info_str *vg_info; /* Vgroup information (for L3 prods)  */
} fmt_str;

/*
 *  This structure is here to carry the l3 bin organization info
 *  # rows, bin size, etc
 */
typedef struct l3_org_str_d {
    int32 numrows; /*  # bin rows */
    int32 bins_eq; /*  # nins at the equator  */
    float64 vsize; /*  degree vertical size of bin  */
    float64 *hsize; /*  degree horizontal size of bin  */
    int32 *start_bin; /*  start bin number for each row  */
    int32 *max_bin; /*  maximun # bins that can be in a row  */
} l3_org_str;
/*
 *  add multipliers for float uncertainty in checking of attributes to
 *  a constant
 */
#define ERR_FRAC_F32 0.0001
#define ERR_FRAC_F64 0.000001
/*
 *  Do some prototype definitions
 */
char *s_parse(char*, int);
int fmt_rd_attr(char *, FILE *, fmt_str *);
int fmt_read(char *, fmt_str *);
char *get_line(char*, int, FILE *, int);
int var_decode(char *, int32, void *, int32, int32);
int get_attr(int32, attr_str, u_data *, int *);
void attr_disp(attr_str, u_data, int);
void chk_str(attr_str, char *, int32);
int fmt_rd_dim(char *, FILE *, fmt_str *);
int fmt_rd_sds(char *, FILE *, fmt_str *);
int fmt_rd_l3vg(char *, FILE *, fmt_str *);
int chk_sds(int32, fmt_str *, int);
int fmt_rd_ras(char *, FILE *, fmt_str *);
int hdf_ras_chk(char *, char *, char *, int, int);
void fmt_exit(int);
void ck_v_l3(int32, fmt_str *);
int32 group_fnd(int32, int32, int32, int32, vg_info_str *);
void chk_sea_grid(int32);
void chk_bin_index(int32, int32 *, int32 *);
void chk_bin_list(int32, int32, int32 *, int32 *);
void chk_l3_prod(int32, int32, char *, char *);
void l3_get_org(int, l3_org_str *);

#ifndef DBL_DIG
#define DBL_DIG 15
#endif

#ifndef FLT_DIG
#define FLT_DIG 6
#endif
