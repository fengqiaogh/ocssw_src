/*
 * prism.h
 *
 *  Created on: June 2015
 *      Author: rhealy
 */

#ifndef SRC_L2GEN_PRISM_H_
#define SRC_L2GEN_PRISM_H_
#include <gsl/gsl_errno.h>
#include <gsl/gsl_spline.h>
#include <gsl/gsl_sort_double.h>
#include <proj.h>

#ifdef  __cplusplus
extern "C" {
#endif

typedef struct prism_l1b_t {
    // info

    int npixels; /**< number of pixels in prism              */
    int nscans; /**< number of scans in prism              */
    int nbands; /**< number of visible bands in prism      */


} prism_l1b_t;

typedef struct prism4ocia_struct {
    int32_t year, day, month, doy, msec, hour, min;
    float sec;
    double stime, etime;
    int32_t npix, nscan, wgs_nscan, wgs_npix;
    float *sena, *senz, *sola, *solz, *utc;
    double *gain, *lon, *lat, scantime;
    double *wave, *fwhm;
    PJ *pj;
    double easting, northing, rotation;
    double pixelSize;
    int utmZone, numBands;
    int interleave, eastbyscan;
    float *Lt;
    float *scale_factor;
    FILE *av_fp;
    gsl_spline *spline;
    gsl_interp_accel *spl_acc;
    float alt;
} prism4ocia_t;


#ifdef  __cplusplus
}
#endif


#endif /* SRC_L2GEN_PRISM_H_ */
