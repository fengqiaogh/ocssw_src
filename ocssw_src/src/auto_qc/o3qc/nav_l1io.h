#ifndef NAVIGATION_H
#define NFLAG_FIX
#define NAVIGATION_H

#include <stdint.h>

#ifndef NAVBLOCK_STRUCT
#define NAVBLOCK_STRUCT

#include "stdint.h"

typedef struct navblockStruct {
    float orb_vec[3];
    float l_vert[3];
    float sun_ref[6];
    float att_ang[3];
    float sen_mat[9]; /*  actually a 3 x 3 array */
    float scan_ell[6];
    int32_t nflag[8];
} navblockType;
#endif /* NAVBLOCK_STRUCT */

#ifndef GEOLOC_STRUCT
#define GEOLOC_STRUCT

typedef struct geolocStruct {
    float *ylat; /* pixel geodetic latitudes		*/
    float *xlon; /* pixel geodetic longitude 		*/
    float *solz; /* pixel solar zenith angle		*/
    float *sola; /* pixel solar azimuth angle		*/
    float *senz; /* pixel sensor zenith angle		*/
    float *sena; /* pixel sensor azimuth angle		*/
} geolocType;
#endif /* GEOLOC_STRUCT */

#endif /* NAVIGATION_H */
