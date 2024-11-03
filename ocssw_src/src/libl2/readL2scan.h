#ifndef READL2SCAN_H
#define READL2SCAN_H

#include <dfutils.h>
#include <stdint.h>

#define byte uint8_t

#define MAXNFILES 10000 /* Increase to 544 09/25/06 JMG - 800 WDR */
#define MAXNUMBERPRODUCTS 2000
#ifdef __cplusplus
extern "C" {
#endif

// set of functions defined in expand 3D module
void l3_l2_conversion(char **inp3, char **out2);
void l3_l2_index(char **inp3, int *start, int *count);
int get_set_flag();
int get_l2_flag_use();
typedef struct l2prod_struct {
    int32_t fileindex;

    char filename[256];
    char oformat[32];

    int32_t nrec;
    int32_t nsamp;

    int16_t syear;
    int16_t sday;
    int32_t smsec;
    int16_t eyear;
    int16_t eday;
    int32_t emsec;
    int32_t orbit;
    char dtype[8];

    int32_t ntilts;
    int16_t tilt_flags[20];
    int16_t tilt_ranges[2][20];

    char *flagnames;
    int32_t flagmask;

    int32_t year;
    int32_t day;
    int32_t msec;

    int32_t *year_cache;
    int32_t *day_cache;
    int32_t *msec_cache;
    double *scantime_cache;

    float *longitude;
    float *latitude;

    float *geonav[6];

    float *lon_cntl;
    float *lat_cntl;
    float *cntl_pnts;
    float *cntl_pnts_cache;
    float *spline_arr;

    int32_t nprod;
    char *prodname[MAXNUMBERPRODUCTS];
    float bv_unscaled[MAXNUMBERPRODUCTS];
    int16_t bv_scaled[MAXNUMBERPRODUCTS];
    int32_t thirdDim[MAXNUMBERPRODUCTS];

    float **l2_data;
    int32_t *l2_flags;
    int32_t *l2_bits;

    byte eng_qual[4];
    byte s_flags[4];
    int32_t nflag[8];

    int32_t geointerp;

    byte *mside;
    byte *detnum;
    int32_t *pixnum;

    // area weighting variables
    float *lat1;     // lat1 for corner of pixel
    float *lon1;     // lon1 for corner of pixel
    float *lat2;     // lat2 for corner of pixel
    float *lon2;     // lon2 for corner of pixel
    byte lat2Valid;  // does lat2 and lon2 contain valid coordinates

    int wavelength[250];
    int wavelengths_3d[250];

    int32_t bandsPerPixel;
    int32_t wavelength_3d;
} l2_prod;

typedef struct meta_l2Struct {
    char *product_name; /* ATTR Product name(file name)         */
    char *title;
    char *data_center;  /* ATTR data_center, processing center  */
    char *mission;      /* ATTR mission                         */
    char *mission_char; /* ATTR Mission Characteristics         */
    char *sensor_name;  /* ATTR sensor name                     */
    char *sensor;       /* ATTR sensor                          */
    char *sensor_char;  /* ATTR instrumentInformation          */
    char *sw_id;        /* ATTR Software ID                     */
    char *infiles;      /* ATTR Input files                     */
    char *stime;        /* ATTR Start time                      */
    char *etime;        /* ATTR End time                        */
    char *ctime;        /* ATTR scene center time               */
    char *ntime;        /* ATTR Node crossing time              */
    char *snode;        /* ATTR Start Node                      */
    char *enode;        /* ATTR End Node                        */
    int orbnum;         /* ATTR orbit number                    */
    char *norad1;       /* ATTR NORAD elements, first line      */
    char *norad2;       /* ATTR NORAD elements, second line     */
    int pix_start;      /* ATTR LAC Pixel Start Number          */
    int pix_sub;        /* ATTR LAC Pixel Subsampling           */
    int ncrec;          /* ATTR scene center scan line          */
    int nfrec;          /* ATTR number of filled scan line      */
    byte ff_mis;        /* ATTR FF missing frames               */
    byte sd_mis;        /* ATTR SDPS missing frames             */
    float flags_pc[32]; /* MFSD % data for each quality flag    */
    char *lat_units;    /* ATTR Latitude units                  */
    char *lon_units;    /* ATTR Longitude units                 */
    float northlat;     /* ATTR Northernmost latitude           */
    float southlat;     /* ATTR Southernmost latitude           */
    float westlon;      /* ATTR Westernmost longitude           */
    float eastlon;      /* ATTR Easternmost longitude           */
    float startclat;    /* ATTR Start Center Latitude           */
    float startclon;    /* ATTR Start Center Longitude          */
    float endclat;      /* ATTR End Center Latitude             */
    float endclon;      /* ATTR End Center Longitude            */
    float nodel;        /* ATTR Orbit node longitude            */
    int ntilts;         /* MFSD Sensor Tilt                     */
    /* Calibration Vgroup */
    short entry_year;
    short entry_day;
    short ref_year;
    short ref_day;
    short ref_minute;
} meta_l2Type;

/* Prototypes */
void free_rowgroup_cache();
void init_rowgroup_cache();
int32_t get_dtype(int32_t dtype, ds_format_t fileformat);
int32_t openL2(const char *fname, const char *plist, l2_prod *l2_str);
int32_t reopenL2(int32_t fileindex, l2_prod *l2_str);
int32_t readL2(l2_prod *l2_str, int32_t ifile, int32_t recnum, int32_t iprod,
               unsigned char *scan_in_rowgroup);
int32_t readlonlat(l2_prod *l2_str, int32_t ifile, int32_t *start, int32_t *edges,
                   unsigned char *scan_in_rowgroup);
int32_t closeL2(l2_prod *l2_str, int32_t ifile);
int32_t freeL2(l2_prod *l2_str);
int32_t findprod(l2_prod *l2_str, char *prodname);
int32_t readL2meta(meta_l2Type *meta_l2, int32_t ifile);
int32_t freeL2meta(meta_l2Type *meta_l2);
int32_t getL3units(l2_prod *l2_str, int32_t ifile, char *l3b_prodname, char *units);

enum L2PixelMode_t { L2PixelOff, L2PixelCorner, L2PixelDelta };

/**
 * Turn on the calculation of the lat/lon deltas.  Used for area weighting.
 *
 * @param val L2PixelOff(0)=off, L2PixelCorner(1)=pixel corners, L2PixelDelta(2)=pixel deltas
 */
void enableL2PixelArea(enum L2PixelMode_t val);

#ifdef __cplusplus
}
#endif
#endif
