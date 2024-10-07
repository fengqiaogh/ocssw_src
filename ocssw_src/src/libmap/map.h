#ifndef MAP_H /* avoid re-inclusion */
#define MAP_H

#include <stdint.h>

#define MAXVAL     255

#define NOMATCH_ERR  -2

/*  Time Conversion Constants  */
#define MSECHOUR 3600000
#define MSECMIN  60000
#define MSECSEC  1000

/* define meta_l3m structure */
typedef struct {
    char *pname;
    char *title;
    char *dcenter;
    char *station;
    float station_lat;
    float station_lon;
    char *mission;
    char *mission_char;
    char *sensor;
    char *sensor_name;
    char *sensor_char;
    char *prodtype;
    char *replace;
    char *softid;
    char *softnm;
    char *softver;
    char *ptime;
    char *infiles;
    char *proccon;
    char *proclog;
    char *flag_use;
    uint8_t eng_q_use[4];
    int16_t bin_syear;
    int16_t bin_sday;
    int16_t bin_eyear;
    int16_t bin_eday;
    char *stime;
    char *etime;
    int16_t syear;
    int16_t sday;
    int32_t smsec;
    int16_t eyear;
    int16_t eday;
    int32_t emsec;
    int32_t orbit;
    int32_t start_orb;
    int32_t end_orb;
    char *mapproj;
    char *lat_units;
    char *lon_units;
    float nlat;
    float slat;
    float elon;
    float wlon;
    float lat_step;
    float lon_step;
    int32_t nbins;
    int32_t nrows;
    int32_t ncols;
    char *parameter;
    char *measure;
    char *units;
    char *scaling_type;
    char *scaling_eqn;
    char *base;
    float slope;
    float intercept;
    float scaled_data_min;
    float scaled_data_max;
    char *scaled_data_type;
    float data_min;
    float data_max;
    int32_t bintype;
} meta_struct;

#endif /* MAP_H */

