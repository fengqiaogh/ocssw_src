/*
 * aviris.h
 *
 *  Created on: May 18, 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_AVIRIS_H_
#define SRC_L2GEN_AVIRIS_H_
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_fit.h>
#include <proj.h>
#define AV_MAXBANDS 224

//typedef struct aviris_l1b_t {
//
//    int     npixels;              /**< number of pixels in AVIRIS        */
//    int     nscans;               /**< number of scans in AVIRIS         */
//    int     nbands;               /**< number of visible bands in AVIRIS */
//
//} aviris_l1b_t;

typedef struct aviris_struct {
    int32_t year, day, month, doy, msec;
    int32_t npix, nscan, wgs_nscan, wgs_npix;
    double *sena, *senz, *sola, *solz, *utc, *lon, *lat;
    float *elev, *alt, lat0, lon0, distmin, distmax;
    double *gain;
    double *wave, *fwhm;
    PJ *pj;
    double easting, northing, rotation;
    double pixelSize;
    int utmZone, numBands;
    int interleave, eastbyscan;
    int have_nav, have_gain;
    char hdrfile[FILENAME_MAX], imgfile[FILENAME_MAX], navfile[FILENAME_MAX], gainfile[FILENAME_MAX];
    FILE *av_fp;
    gsl_spline *spline;
    gsl_interp_accel *spl_acc;
    int isnetcdf;
} aviris_t;

#endif /* SRC_L2GEN_AVIRIS_H_ */
