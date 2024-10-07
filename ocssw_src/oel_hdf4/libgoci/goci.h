/**
 *  @file libgoci/goci.h
 *  @brief GOCI file format reader
 *  @author Paul Martinolich
 *  @author Naval Research Laboratory
 */

#ifndef GOCI_H
#define GOCI_H

#include <proj.h>

typedef struct goci_l1b_t {
    // info

    int npixels; /**< number of pixels in GOCI L1B             */
    int nscans; /**< number of scans in GOCI L1B              */
    int nbands; /**< number of visible bands in GOCI L1B      */

    PJ *pj;
    double startX;
    double startY;
    double deltaX;
    double deltaY;
    float sat_pos[3]; /**< position of satellite */

    // HDF5 values

    hid_t fid; /**< HDF5 file descriptor                     */
    hid_t gid; /**< HDF5 "/" root group descriptor           */
    hid_t dset_id[8]; /**< Dataset IDs for all 8 bands              */
    int32_t slot_nav_avail; /**< is GOCI slot navigation available? */
    float *slot_rel_time; /**< slot mean time relative to start */
    unsigned char *slot_asg; /**< pixel by pixel slot assignment */

} goci_l1b_t;

/*
 *  GOCI slot nav table from L1B
 */
typedef struct slot_nav_str_def {
    int32_t band_num;
    int32_t slot_num;
    float rel_time;
    float sc_att[3];
    float xo;
    float yo;
    float xs;
    float ys;
    float xpo;
    float ypo;
    float xps;
    float yps;
    int32_t num_a_parm;
    float a_parm[16];
    int32_t num_b_parm;
    float b_parm[16];
    int32_t num_c_parm;
    float c_parm[16];
    int32_t num_d_parm;
    float d_parm[16];
    int32_t num_ap_parm;
    float ap_parm[16];
    int32_t num_bp_parm;
    float bp_parm[16];
    int32_t num_cp_parm;
    float cp_parm[16];
    int32_t num_dp_parm;
    float dp_parm[16];
} slot_nav_str;

int goci_proj4_open(goci_l1b_t *l1b);
void goci_proj4_convert(goci_l1b_t *l1b, int numPoints, double *x, double *y);
void goci_proj4_close(goci_l1b_t *l1b);

int goci_l1b_open(const char *src_path, goci_l1b_t **goci_l1b);
int goci_l1b_close(goci_l1b_t *goci_l1b);
int goci_l1b_get_date(goci_l1b_t *goci_l1b, char *tim_str, int *year,
        int *month, int *day);
int goci_l1b_get_time(goci_l1b_t *goci_l1b, char *tim_str, int *hour,
        int *min, int *sec);
int goci_l1b_get_band(goci_l1b_t *goci_l1b, int band, int line,
        uint32_t *buf);
int32_t goci_slot_init(hid_t file_id, hsize_t *dims, float *slot_rel_time,
        unsigned char *slot_asg, int32_t *slot_nav_avail);
int32_t goci_slot_nav(int32_t ipix, int32_t ilin, int32_t bnd, int32_t itile,
        slot_nav_str *slot_nav, int32_t nbnd, int32_t nslot,
        unsigned char *bnd_tile_lut, float *nradsq);
unsigned char goci_slot_time(int32_t ipix, int32_t ilin, goci_l1b_t *goci_l1b,
        float *rel_sec);

#endif
