#ifndef ANC_H /* avoid re-inclusion */
#define ANC_H

#include <stdint.h>
#include <stdio.h>
#include <math.h>
#include <string.h>
#include <ctype.h>

#define LAT  721
#define LON  1441

#define MAXVAL 255
#define null '\0'

#define WIND_U    0
#define WIND_V    1
#define PRESSURE  2
#define HUMIDITY  3
#define OZONE     4
#define TAUAER    6
#define ANGSTROM  7
#define TAUCLD    8
#define CFCLD     9

#define CORNERS   4
#define NPARMS    6
#define NFILE     2
#define F1        0
#define F2        1

/*  Time Conversion Constants  */
#define MSECHOUR 3600000
#define MSECMIN  60000
#define MSECSEC  1000

#define BS_INCR 0       /* increasing order */
#define BS_DECR 1       /* decreasing order */

#define SPATIAL 110
#define SPATIAL_TEMPORAL 112

#define GEOMETRY_VDATA   "Equal-Angle SDS"

#define WIND_U_LABEL   "Zonal Wind"
#define WIND_V_LABEL   "Meridional Wind"
#define PRESSURE_LABEL   "Atmospheric Pressure"
#define HUMIDITY_LABEL   "Relative Humidity"
#define OZONE_LABEL   "Total Ozone"

#define CLIM_DATA                 "Average"

#define QC_LABEL   "QC flag"
#define PARM_LABEL              "Observations"

#define SYEAR  "Start Year"
#define SDAY  "Start Day"
#define SMSEC  "Start Millisec"
#define EYEAR  "End Year"
#define EDAY  "End Day"
#define EMSEC  "End Millisec"

#define VSIZE    "Latitude Step"
#define HSIZE   "Longitude Step"
#define MAXNORTH  "Northernmost Latitude"
#define MAXSOUTH "Southernmost Latitude"
#define MAXWEST  "Westernmost Longitude"
#define MAXEAST  "Easternmost Longitude"
#define SWPT_LAT "SW Point Latitude"
#define SWPT_LON "SW Point Longitude"

extern char ERR_MSG[1024];

#ifndef DATANAMES
#define DATANAMES
static const char * const data_sdsnames[NPARMS] = {
    "z_wind",
    "m_wind",
    "press",
    "p_water",
    "ozone",
    "rel_hum"
};

static const char * const QC_sdsnames[NPARMS] = {
    "z_wind_QC",
    "m_wind_QC",
    "press_QC",
    "p_water_QC",
    "ozone_QC",
    "rel_hum_QC"
};

static const char * const clim_vgps[12] = {
    "January",
    "February",
    "March",
    "April",
    "May",
    "June",
    "July",
    "August",
    "September",
    "October",
    "November",
    "December"
};

static const char * const clim_datasets[NPARMS] = {
    "z_wind_mean",
    "m_wind_mean",
    "press_mean",
    "p_water_mean",
    "ozone_mean",
    "rel_hum_mean"
};
#endif /* DATANAMES */
/*
#ifndef FILESTRUCT
#define FILESTRUCT

#define FLIMIT  12
struct file_time {
        char     *fn;
        float64  sdate;
        float64  stime;
        float64  edate;
        float64  etime;
 } ftime[FLIMIT];
#endif */ /* FILESTRUCT */


/**** fortran prototypes */

void dataintp_(float* in_latlon, float* lat_list, float* lon_list,
        float* data_list1, double* DT1, float* data_list2, double* DT2,
        int32_t* ipt, int32_t* nband, float(*) [2], float* def,
        int32_t* intporder, float* dummy, float* dataout, float* unc,
        int32_t* int_bad, int32_t* row, int32_t* col);


int dataintp(float in_latlon[2], float* lat_list, float* lon_list,
               float* data_list1, double* DT1, float* data_list2, double* DT2,
               int32_t* ipt, int32_t* nband, float rng[2], float* def,
               int32_t* intporder, float* dummy, float* dataout, float* unc,
               int32_t* int_bad, int32_t* row, int32_t* col);

void julian_(double* dtin, double* d_jd);



/**** internal c prototypes */

int32_t rdancattr(int32_t sdfid, char *attr_name, void *buf);



#endif /* ANC_H */
