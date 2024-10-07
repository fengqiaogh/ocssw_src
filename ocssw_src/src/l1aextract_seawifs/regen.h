#ifndef REGEN_H
#define REGEN_H

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <time.h>
#include "hdf.h"
#include "mfhdf.h"

#define NBANDS 8
#define MSECHOUR 3600000
#define MSECMIN  60000
#define MSECSEC  1000
#define PFLAG   16

#define MAXATTRS        255
#define MAXVAL  255
#define MALLOC_ERR   -5
#define MAXVGPS    6 

#ifndef ERR_MSG
extern char ERR_MSG[255];
#endif

#define  DTYPE          "Data Type"
#define  NSAMP          "Pixels per Scan Line"
#define  NREC           "Number of Scan Lines"
#define  TILTFLAGS  "tilt_flags"
#define  TILTRANGES     "tilt_ranges"
#define  TILTLATS   "tilt_lats"
#define  TILTLONS   "tilt_lons"
#define  PTIME  "Processing Time"
#define  SYEAR          "Start Year"
#define  SDAY           "Start Day"
#define  SMSEC          "Start Millisec"
#define  EYEAR          "End Year"
#define  EDAY           "End Day"
#define  EMSEC          "End Millisec"
#define  PIX_START      "LAC Pixel Start Number"
#define  PIX_SUB        "LAC Pixel Subsampling"
#define  CTIME          "Scene Center Time"
#define  SCLAT          "Scene Center Latitude"
#define  SCLON          "Scene Center Longitude"
#define  SCSOLZ         "Scene Center Solar Zenith"
#define  ULLAT          "Upper Left Latitude"
#define  ULLON          "Upper Left Longitude"
#define  URLAT          "Upper Right Latitude"
#define  URLON          "Upper Right Longitude"
#define  LLLAT          "Lower Left Latitude"
#define  LLLON          "Lower Left Longitude"
#define  LRLAT          "Lower Right Latitude"
#define  LRLON          "Lower Right Longitude"
#define  NORTHLAT       "Northernmost Latitude"
#define  SOUTHLAT       "Southernmost Latitude"
#define  WESTLON        "Westernmost Longitude"
#define  EASTLON        "Easternmost Longitude"
#define  BCLAT  "Start Center Latitude"
#define  ECLAT  "End Center Latitude"
#define  BCLON  "Start Center Longitude"
#define  ECLON  "End Center Longitude"
#define  NCREC          "Scene Center Scan Line"
#define  NFREC          "Filled Scan Lines"
#define  FF_MIS         "FF Missing Frames"
#define  SD_MIS         "SDPS Missing Frames"
#define  PCT_FLAG       "Flag Percentages"

#define  L1SATG1        "Gain 1 Saturated Pixels"
#define  L1SATG2        "Gain 2 Saturated Pixels"
#define  L1NSATG1       "Gain 1 Non-Saturated Pixels"
#define  L1NSATG2       "Gain 2 Non-Saturated Pixels"
#define  L1ZEROS        "Zero Pixels"
#define  L1MEANR1       "Mean Gain 1 Radiance"
#define  L1MEANR2       "Mean Gain 2 Radiance"

#ifndef TILT_STRUCT
#define TILT_STRUCT

typedef struct tilt_struct {
    int ntilts; /* MFSD Number of scene tilt states*/
    short tilt_flags[20]; /* MFSD Tilt indicators         */
    short tilt_ranges[20][2]; /* MFSD ranges of scene tilt states */
    float tilt_lats[20][2][2]; /* MFSD Lats of tilt-range pts  */
    float tilt_lons[20][2][2]; /* MFSD Lons of tilt-range pts  */
} tilt_Type;
#endif /* TILT_STRUCT */


#ifndef NAV_STRUCT
#define NAV_STRUCT

typedef struct nav_struct {
    float *orb_vec; /* MFSD Orbit positionn vector  */
    float *l_vert; /* Local vertical vector in ECEF*/
    float *sun_ref; /* Reference Sun vector in ECEF */
    float *att_ang; /* Computed yaw, roll, pitch    */
    float *sen_mat; /* ECEF-to-sensor-frame matrix  */
    float *scan_ell; /* Scan-track ellipse coeffs    */
    int32 *nflag; /* Navigation flags             */
} NavType;
#endif /* MET_STRUCT */

#ifndef GEO_STRUCT
#define GEO_STRUCT

typedef struct geonav_struct {
    float *slat; /* MFSD Scan start pixel latitude*/
    float *slon; /* MFSD Scan start pixel longitude*/
    float *clat; /* MFSD Scan center pixel lat   */
    float *clon; /* MFSD Scan center pixel long  */
    float *elat; /* MFSD Scan end pixel latitude */
    float *elon; /* MFSD Scan end pixel longitude*/
    float *csol_z; /* MFSD Scene Center Solar Zenith*/
} GeoType;
#endif /* GEO_STRUCT */

#ifndef FM_STRUCT
#define FM_STRUCT

typedef struct fm_struct {
    int32 satg1[8];
    int32 satg2[8];
    int32 nsatg1[8];
    int32 nsatg2[8];
    int32 zeroes[8];
    float meanr1[8];
    float meanr2[8];
} FilemetricsType;
#endif /* FM_STRUCT */

#endif /* REGEN_H */
