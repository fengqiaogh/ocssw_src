#ifndef _L1_INPUT_H
#define _L1_INPUT_H

#include <stdint.h>
#include <stdio.h>

#define L1_PRODSTRLEN 2048
#define L1_MAXPROD 1000
#define L1_NFLAGS 32

#define NBANDSIR       8

#ifdef __cplusplus
extern "C" {
#endif


typedef struct l1_input_struct {

    // file reading control params
    char calfile[FILENAME_MAX];
    char xcal_file[FILENAME_MAX];
    char btfile[FILENAME_MAX];
    char cld_msk_file[FILENAME_MAX];
    char viirscalparfile[FILENAME_MAX];
    char pversion[1024];
    char input_parms[64000];  // lim is 65K, just do this
    char input_files[6144];

    int32_t rad_opt; /* radcor switch for MERIS smile correction */
    int32_t geom_per_band; /* 0 - use nominal geometry sen/sol_a/a
                                 1 - use band-specific values for instruments
                                 that have band specific geometry */
    int32_t xcal_nwave; /* number of wavelengths to which xcal applied */
    int32_t *xcal_opt; /* xcal option per band              */
    float *xcal_wave; /* sensor wavelengths to which xcal applied */
    int32_t resolution; /* process at this nadir pixel res (meters) */
                        /* 250, 500, 1000, -1=native (modis only)   */
    int32_t newavhrrcal; /* new avhrr calibration equation */
    int32_t sl_pixl; /* seawifs straylight pixel limit           */
    float sl_frac; /* seawifs straylight Ltyp fraction         */
    float ch22detcor[10]; /* channel 22 detector corrections  */
    float ch23detcor[10]; /* channel 23 detector corrections  */
    float cirrus_thresh[2]; /* cirrus reflectance thresholds    */
    float albedo; /* cloud reflectance threshold      */
    float cloud_wave; /* cloud test wavelength            */
    float cloud_eps; /* cloud reflectance ratio          */
    float glint; /* glint threshold                  */
    float extreme_glint;/* extreme glint threshold   */
    float sunzen; /* solar zenith angle threshold     */
    float satzen; /* sensor zenith angle threshold    */
    float hipol; /* high polarization threshold      */
    float *gain; /* Vicarious calibration gain       */
    float *offset; /* Vicarious calibration offset     */

    int32_t spixl; /* starting pixel no. of the input (1-rel)  */
    int32_t epixl; /* ending pixel no. of the input (1-rel)    */
    int32_t dpixl; /* pixel subsampling increment              */
    int32_t sline; /* starting line no. of the input (1-rel)   */
    int32_t eline; /* ending line no. of the input (1-rel)     */
    int32_t dline; /* line subsampling increment               */

    int32_t outband_opt; /* 1=apply seawifs out-of-band correction   */
    int32_t evalmask;
    int32_t landmask; /* 0=off, 1=on */
    int32_t bathmask; /* 0=off, 1=on */
    int32_t cloudmask; /* 0=off, 1=on */
    int32_t glintmask; /* 0=off, 1=on */
    int32_t sunzenmask; /* 0=off, 1=on */
    int32_t satzenmask; /* 0=off, 1=on */
    int32_t hiltmask; /* 0=off, 1=on */
    int32_t stlightmask; /* 0=off, 1=on */
    int32_t cloud_mask_opt; /* 0=cloud_flag, 1=cloud_flag_dilated */
    
} l1_input_t;



float* calloc_nbandsf(int32_t nbands, float *nbarray, float init_val);
int32_t* calloc_nbandsi32t(int32_t nbands, int32_t *nbarray, int32_t init_val);
int* calloc_nbandsi(int32_t nbands, int *nbarray, int init_val);



#ifdef __cplusplus
}
#endif

extern l1_input_t *l1_input; /* L1 input parameters structure         */


#endif
