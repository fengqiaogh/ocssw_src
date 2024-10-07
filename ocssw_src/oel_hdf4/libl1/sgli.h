/*
 * sgli.h
 *
 *  W. Robinson, SAIC, 3 Oct 2016
 *
 */

#ifndef SRC_SGLI_SGLI_H_
#define SRC_SGLI_SGLI_H_

#include <gsl/gsl_errno.h>
#include <gsl/gsl_interp2d.h>
#include <gsl/gsl_spline2d.h>

#include "h5io.h"
#include "filehandle.h"

typedef struct sgli_struct {
    char swir_ir_file[FILENAME_MAX];
} sgli_t;

/* This is a grid-specific definition - 2 are expected for SGLI */
typedef struct grid_res_str_def {
    double *xa, *ya;  /* point to 1 of 2 grid point arrays for, x, y of
                         this grid  */
    int32_t npix_tie, nscn_sub_tie, tie_st_lin;  /* tie point array
                       x, y size and actual start in full array */
    int32_t resamp;  /* resampling to the dominent resolution */
} grid_res_str;

/* This is a sensor band-dependent geometry structure, per-band  */
typedef struct band_geom_str_def {
    h5io_str dsid[2]; /* for sensor as solar has no band dependence */
    float scale[2];
    float offset[2];
    grid_res_str *grd_desc;  /* specific grid resolution structure for this one */
    gsl_spline2d *int_id_sen[3];
    char *qual;
} band_geom_str;

#endif /* SRC_SGLI_SGLI_H_ */
