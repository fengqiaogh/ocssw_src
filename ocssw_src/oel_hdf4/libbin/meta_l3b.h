#ifndef METAL3B_H /* avoid re-inclusion */
#define METAL3B_H

#include <dfutils.h>
#include <hdf5.h>

#define SM_ATTRSZ 1024
#define MD_ATTRSZ 10000
#define LG_ATTRSZ 65535        /* MAX_ORDER is defined in HDF as 65535 */

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct {
    char product_name[SM_ATTRSZ];
    char title[SM_ATTRSZ];
    int sensorID;
    char sensor_name[SM_ATTRSZ];
    char data_center[SM_ATTRSZ];
    char mission[SM_ATTRSZ];
    char mission_char[SM_ATTRSZ];
    char sensor[SM_ATTRSZ];
    char sensor_char[SM_ATTRSZ];
    char station[SM_ATTRSZ];
    float station_lat;
    float station_lon;
    char units[MD_ATTRSZ];
    char prod_type[SM_ATTRSZ];
    char pversion[SM_ATTRSZ];
    char replace[SM_ATTRSZ];
    char soft_name[SM_ATTRSZ];
    char soft_ver[SM_ATTRSZ];
    char ptime[SM_ATTRSZ];
    char proc_con[MD_ATTRSZ];
    char input_parms[LG_ATTRSZ];
    char flag_names[SM_ATTRSZ];
    char infiles[LG_ATTRSZ];
    double startTime;
    double endTime;
    int32_t orbit;
    int32_t start_orb;
    int32_t end_orb;
    char lat_units[SM_ATTRSZ];
    char lon_units[SM_ATTRSZ];
    int64_t data_bins;
    float pct_databins;
    char binning_scheme[SM_ATTRSZ];
    float north;
    float south;
    float east;
    float west;
    double resolution; // in meters
    char doi[SM_ATTRSZ];
    char keywords[SM_ATTRSZ];
} meta_l3bType;

void calculate_temporal_range(meta_l3bType *meta_l3b);
int write_l3b_meta_hdf4(int32_t sd_id, meta_l3bType *meta_l3b);
int write_l3b_meta_netcdf4(idDS ds_id, meta_l3bType *meta_l3b, int write64bit);

int32_t rdattr(int32_t sdfid, char *attr_name, void *buf);
int32_t getattrsz(int32_t id, char *attr_name, int32_t *nt, int32_t *count);
int32_t read_l3b_meta_hdf4(int32_t sdfid, meta_l3bType *meta_l3b);
int32_t read_l3b_meta_hdf5(hid_t grp0, meta_l3bType *meta_l3b);
int32_t read_l3b_meta_netcdf4(int ncid, meta_l3bType *meta_l3b);

#ifdef  __cplusplus
}
#endif

#endif /* METAL3B_H */
