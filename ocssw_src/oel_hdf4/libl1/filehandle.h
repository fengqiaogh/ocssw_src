#ifndef _FILEHANDLE_H
#define _FILEHANDLE_H

#include <stdio.h>
#include <string.h>
#include <filetype.h>
#include <sensorInfo.h>
#include <productInfo.h>
#include "l2prod_struc.h"

#define FORWARD      0
#define INVERSE_ZERO 1
#define INVERSE_NLW  2
#define INVERSE_LW   3

#define READ  0
#define WRITE 1

#define L1_PRODSTRLEN 2048
#define L1_MAXPROD 1000
#define L1_NFLAGS 32

#define NBANDSIR       8

#ifdef __cplusplus
extern "C" {
#endif

typedef struct filehandle_struct {
    char name[FILENAME_MAX];
    file_type format;
    int32_t sensorID;
    int32_t subsensorID;
    char spatialResolution[10];
    int32_t length;
    int32_t spix; /* start pixel (0-based)                  */
    int32_t epix; /* end pixel (0-based)                    */
    int32_t npix;
    int32_t nscan;
    int32_t nbands;
    int32_t nbandsir;
    int32_t nlvl;
    int32_t n_refl_loc;
    int32_t n_cloud_phase;
    int32_t *bindx; /* index to closest seawifs band          */
    int32_t ndets;
    int32_t mode;
    char l2prod[L1_PRODSTRLEN]; /* list of L2 products to be included     */
    char def_l2prod[L1_PRODSTRLEN]; /* list of default L2 products      */
    int32_t sd_id; /* hdf file id for the opened output file */
    int32_t tot_prod; /* total # of L2 products to be created   */
    char l2_prod_names[L1_MAXPROD][32];
    l2prodstr* prodptr; /* array of product structures */
    productInfo_t** productInfos; // pointers to product info structures

    char *geofile;
    char *gmpfile;
    float orbit_node_lon;
    int32_t orbit_number;
    double node_crossing_time;
    int32_t flag_cnt[L1_NFLAGS];
    int32_t terrain_corrected;
    int32_t sv_with_moon;
    int32_t grp_id[8]; // netCDF group IDs

    int32_t *iwave;
    float *fwave;
    float *fwhm;
    float *Fobar;
    float *Fonom;
    float *Tau_r;
    float *k_oz;
    float *k_no2;
    float *aw;
    float *bbw;

    void *private_data;

} filehandle;

float* calloc_nbandsf(int32_t nbands, float *nbarray, float init_val);
int32_t* calloc_nbandsi32t(int32_t nbands, int32_t *nbarray, int32_t init_val);
int* calloc_nbandsi(int32_t nbands, int *nbarray, int init_val);

void filehandle_init(filehandle *file);



#ifdef __cplusplus
}
#endif

#endif
