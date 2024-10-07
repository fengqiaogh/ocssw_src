#ifndef L2BRSGEN_H /* avoid re-inclusion */
#define L2BRSGEN_H

#include  <stdio.h>
#include  <string.h>
#include        <time.h>

#include        "hdf.h"
#include        "mfhdf.h" 

#include <clo.h>
#include <readL2scan.h>
#include <get_product_table.h>

//#include        "hdfhdr.h"
//#include        "hdfmac.h"
//#include        "usrhdr.h"
#include        "usrmac.h"

#define VERSION "2.0"

#define IOERR           -1
#define RDERR           -2
#define MEMERR          -3

#define MAXVAL  255

#define CHL_START       (l2_data + 8*nsamp)
#define CHL_END         (l2_data + 9*nsamp)
#define LAT_START       (l2_str.latitude)
#define LAT_END         (l2_str.latitude + max_samp_used+1)
#define LON_START       (l2_str.longitude)
#define LON_END         (l2_str.longitude + max_samp_used+1)

#define MAXTILTS 20
#define NFLAGS   16

/*  Time Conversion Constants  */
#define MSECHOUR 3600000
#define MSECMIN  60000
#define MSECSEC  1000

#define LONGNAME "long_name"
#define RANGE  "valid_range"

/*  Global (file-level) Attributes  */
/*  Constants defined for attributes defined in sec.5.3.1 of product specs */

#define L2B_PNAME           "Product Name"
#define L2B_TITLE     "Title"
#define LEGEND  "Legend"
#define TITLE_VAL   " Level-2 Browse Data"
#define DCENTER     "Data Center"
#define DCENTER_VAL  "NASA/GSFC SeaWiFS Data Processing Center"
#define L2BRS_MISSION     "Mission"
#define MSNCHAR         "Mission Characteristics"
#define SENSOR_NAME     "Sensor Name"
#define SENSOR     "Sensor"
#define SNSCHAR  "Sensor Characteristics"


#define REPLACES    "Replacement Flag"
#define SOFTID      "Software ID"
#define SOFTID_VAL  "?????"
#define SOFT_NAME       "Software Name"
#define SOFT_VER        "Software Version"
#define PTIME       "Processing Time"
#define L2BRS_INFILES   "Input Files"
#define PROC_CON    "Processing Control"
#define PROC_LOG    "Processing Log"

/*  Constants defined for attributes defined in sec.5.3.2 of product specs */

#define PINFILES   "Parent Input Files"
#define PPNAME     "Parent Product Name"
#define DTYPE      "Data Type"
#define NSAMP      "Parent Pixels per Scan Line"
#define NREC       "Parent Number of Scan Lines"
#define SNCNTR    "Scene Center Scan Line"
#define NFREC      "Filled Scan Lines"
#define FFLAG1    "FF Missing Frames"
#define FFLAG2     "SDPS Missing Frames"
#define L2BRS_PCTFLAG    "Flag Percentages"

/*  Constants defined for attributes defined in sec.5.3.3 of product specs */
#define STIME       "Start Time"
#define END_TIME       "End Time"
#define CTIME       "Scene Center Time"
#define NTIME     "Node Crossing Time"
#define SYEAR       "Start Year"
#define SDAY        "Start Day"
#define SMSEC       "Start Millisec"
#define EYEAR       "End Year"
#define EDAY        "End Day"
#define EMSEC       "End Millisec"
#define SNODE     "Start Node"
#define ENODE     "End Node"
#define ORBNUM      "Orbit Number"
#define NORAD1      "NORAD Line 1"
#define NORAD2      "NORAD Line 2"

/*  Constants defined for attributes defined in sec.5.3.4 of product specs */
#define LATUNITS      "Latitude Units"
#define LATUNITS_VAL   "degrees North"
#define LONUNITS      "Longitude Units"
#define LONUNITS_VAL   "degrees East"
#define CLAT          "Scene Center Latitude"
#define CLON          "Scene Center Longitude"
#define SCSOL_Z   "Scene Center Solar Zenith"
#define ULLAT   "Upper Left Latitude"
#define ULLON   "Upper Left Longitude"
#define URLAT   "Upper Right Latitude"
#define URLON   "Upper Right Longitude"
#define LLLAT   "Lower Left Latitude"
#define LLLON   "Lower Left Longitude"
#define LRLAT   "Lower Right Latitude"
#define LRLON   "Lower Right Longitude"
#define NLAT          "Northernmost Latitude"
#define SLAT       "Southernmost Latitude"
#define WLON       "Westernmost Longitude"
#define ELON       "Easternmost Longitude"
#define STCLAT   "Start Center Latitude"
#define STCLON   "Start Center Longitude"
#define ENDCLAT   "End Center Latitude"
#define ENDCLON   "End Center Longitude"
#define NODEL         "Orbit Node Longitude"

/*  Constants defined for attributes defined in sec.5.3.4 of product specs */
#define PARAM        "Parameter"
#define PARAM_VAL    "Chlorophyll a concentration"
#define UNITS        "Units"
#define UNITS_VAL    "mg m^-3"
#define PX_START     "Start Pixel"
#define PX_END     "End Pixel"
#define LAC_PX_ST    "LAC Pixel Start Number"
#define PX_SUBSAMP   "Pixel Subsampling Rate"
#define LAC_PX_SUBSAMP   "LAC Pixel Subsampling"
#define PX_NUM       "Pixels per Scan Line"
#define SC_START     "Start Scan"
#define SC_END     "End Scan"
#define SC_SUBSAMP   "Scan Subsampling Rate"
#define SC_NUM       "Number of Scan Lines"
#define PX_LL_NUM    "Pixel Coordinates"
#define SC_LL_NUM    "Scan Coordinates"
#define SC_TYPE      "Scaling"
#define SC_TY_VAL    "logarithmic"
#define SC_EQN       "Scaling Equation"
#define SC_EQN_VAL  "Base**((Slope*brs_data) + Intercept) = chlorophyll a"
#define BASE      "Base"
#define BASE_VAL     10.0
#define SLOPE        "Slope"
#define SLOPE_VAL    0.015
#define INTERCEPT   "Intercept"
#define INTERCEPT_VAL    -2.0
#define SCALE_OFF    "Scale Offset"

/*  Constants defined for attributes defined in sec.5.4.1 of product specs */
#define PX_LL_FIRST   "px_ll_first"
#define PX_LL_FST_ATTR   "Lat/lon of pixels along first scan line"
#define PX_LL_LAST    "px_ll_last"
#define PX_LL_LST_ATTR   "Lat/lon of pixels along last scan line"
#define SC_LL_FIRST   "sc_ll_first"
#define SC_LL_FST_ATTR   "Lat/lon of starts of scan lines"
#define SC_LL_LAST    "sc_ll_last"
#define SC_LL_LST_ATTR   "Lat/lon of ends of scan lines"

/*  Constants defined for attributes defined in sec.5.4.2 of product specs */
#define NTILTS   "ntilts"
#define NTILTS_NAME  "Number of scene tilt states"
#define T_FLAGS       "tilt_flags"
#define T_FLAGS_NAME   "Tilt indicators"
#define T_RANGES      "tilt_ranges"
#define T_RANGES_NAME  "Scan-line number ranges of scene tilt states"
#define T_LATS        "tilt_lats"
#define T_LATS_NAME  "Latitudes of tilt-range scan line end points"
#define T_LONS        "tilt_lons"
#define T_LONS_NAME  "Longitudes of tilt-range scan line end points"

/* Constants defined for navigation data */
#define ORBVEC           "orb_vec"
#define ORBVEC_NAME      "Orbit position vector at scan line time"
#define ORBVEC_UNITS     "kilometers"
#define LVERT            "l_vert"
#define LVERT_NAME       "Local vertical vector in ECEF frame"
#define SUNREF           "sun_ref"
#define SUNREF_NAME      "Reference Sun vector in ECEF frame"
#define ATTANG           "att_ang"
#define ATTANG_NAME      "Computed yaw, roll, pitch"
#define SENMAT           "sen_mat"
#define SENMAT_NAME      "ECEF-to-sensor-frame matrix"
#define SCANELL          "scan_ell"
#define SCANELL_NAME     "Scan-track ellipse coefficients"
#define NFLAG            "nflag"
#define NFLAG_NAME       "Navigation flags"
#define CNTLPTCOLS       "cntl_pt_cols"
#define CNTLPTCOLS_NAME  "Control point column values"

#define CNTLPTLAT        "latitude"
#define CNTLPTLAT_NAME   "Latitude values"
#define CNTLPTLON        "longitude"
#define CNTLPTLON_NAME   "Longitude values"


extern char ERR_MSG[1024];


int l2brsgen_init_options(clo_optionList_t* list);
int l2brsgen_read_options(clo_optionList_t* list, int argc, char* argv[],
        l2_prod *l2_str, meta_l2Type *meta_l2);
int32 put_l2brs(char *l2brs_path, char *replaces, char *ptime, char *infiles,
        int32 px_start, int32 px_end, int32 px_subsamp, int32 brs_nsamp,
        int32 sc_start, int32 sc_end, int32 sc_subsamp, int32 brs_nrec,
        char *l2brs_name, float32 *l2brs_data, int32 *l2brs_flags,
        char *flag_names, int32_t * l2_bits, char *mskflg, unsigned char *palette,
        float32 *px_ll_first, float32 *px_ll_last, float32 *sc_ll_first,
        float32 *sc_ll_last, char *proc_con, int16 syear, int16 sday,
        int32 smsec, int16 eyear, int16 eday, int32 emsec, char *dtype,
        int32 nrec, int32 nsamp, int32 ntilts, short *tilt_flags,
        int16 *tilt_ranges, int16 *cntl_pt_lat, int16 *cntl_pt_lon,
        meta_l2Type *meta_l2, product_table_t *ptable_rec, const char* oformat, int32 apply_pal);

void write_attrs(int32 sdfid, char *l2brs_path, char *replaces, char *ptime,
        char *infiles, int32 px_start, int32 px_end, int32 px_subsamp,
        int32 brs_nsamp, int32 sc_start, int32 sc_end, int32 sc_subsamp,
        int32 brs_nrec, char *l2brs_name, char *proc_con, int16 syear,
        int16 sday, int32 smsec, int16 eyear, int16 eday, int32 emsec,
        char *dtype, int32 nrec, int32 nsamp, meta_l2Type *meta_l2);

int32 write_image(char *l2brs_path, unsigned char *l2brs_data, int32 brs_nsamp,
        int32 brs_nrec, uint8 *palette);

int32 write_SDS(int32 sdfid, char *label, int32 ntype, int32 rank,
        int32 *dimsizes, int32 *start, void *buf);

int32 write_tilt_sets(int32 fid, int32 sdfid, int32 ntilts, int16 *tilt_flags,
        int16 *tilt_ranges);

int32 write_nav_sets(int32 fid, int32 sdfid, int32 brs_nrec, int32 brs_nsamp,
        int16 *cntl_pt_lat, int16 *cntl_pt_lon);


#endif /* L2BRSGEN_H */
