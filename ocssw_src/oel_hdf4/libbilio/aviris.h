/*
 * aviris.h
 *
 *  Created on: May 18, 2015
 *      Author: rhealy
 */

#ifndef BILIO_AVIRIS_H_
#define BILIO_AVIRIS_H_
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>
#include <gsl/gsl_fit.h>
#include <proj.h>

#ifdef __cplusplus
extern "C" {
#endif

//typedef struct aviris_l1b_t {
//
//int     npixels;              /**< number of pixels in AVIRIS        */
//int     nscans;               /**< number of scans in AVIRIS         */
//int     nbands;               /**< number of visible bands in AVIRIS */

//} aviris_l1b_t;

typedef struct aviris4orca_struct {
    int32_t year, day, month, doy, msec, hour, min;
    float sec;
    int32_t npix, nscans, wgs_nscan, wgs_npix;
    //    double **sena, **senz, **sola, **solz, **utc;
    float *sena, *senz, *sola, *solz, *utc;
    //    double *lat, *lon, *elev;
    float *elev, lat0, lon0, distmin, distmax;
    float *Lt;
    float *scale_factor, *alt;
    float *wave, *fwhm;
    PJ *pj;
    float easting, northing, rotation;
    float pixelSize;
    double *gain, *lon, *lat, scantime;
    int utmZone, numBands;
    int interleave, eastbyscan;
    int have_nav, have_gain;
    FILE *av_fp;
    gsl_spline *spline;
    gsl_interp_accel *spl_acc;
} aviris4ocia_t;



int close_aviris(aviris4ocia_t *data);
aviris4ocia_t* open_aviris(char *filename, char *imgfile, char *navfile, char *gainfile, aviris4ocia_t **data);
int read_aviris(aviris4ocia_t *data, int32_t recnum);
int checkAvProcessFile(char *filename, char *hdrfile, char *imgfile, char *navfile, char *gainfile, int itemsize);


#ifdef __cplusplus
}
#endif

#endif /* BILIO_AVIRIS_H_ */
